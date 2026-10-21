# ==============================================================================
# run_model.R — Development/exploration script (2023-2024)
#
# Runs a single baseline simulation and produces diagnostic plots.
# For paper figures, use run_figures.R instead.
#
# Pipeline:
#   1. Source config, data prep, model core (R/00-02)
#   2. Burn-in (40yr cycled climate) -> endemic state
#   3. Forward experiment (2023-2024 actual climate)
#   4. Diagnostic plots: climate, FOI, time series, spatial maps
#   5. Optional hindcast block (2000-2024) if EXPERIMENT_MODE == "hindcast"
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(Matrix)
  library(tibble)
  library(zoo)
  library(patchwork)
  library(stringr)
  library(scales)
  library(igraph)
  library(viridis)
  library(readr)
})

# ==============================================================================
# Source components
# ==============================================================================

# EXPERIMENT_MODE controls whether the hindcast block runs.
# Set to "hindcast" to run 2000-2024 validation, or "contemporary" for 2023-2024 only.
# Must be set BEFORE sourcing 01_data_prep.R so hindcast matrices are built only when needed.
EXPERIMENT_MODE <- "contemporary"

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_model.R")
source("R/03_helpers.R")

message("\n--- Configuration ---\n")
cat("  Experiment mode: ", EXPERIMENT_MODE, "\n")
cat("  Temperature dependence: ", USE_TEMP_DEPENDENCE, "\n")
cat("  Market movement: ", USE_MOVEMENT, "\n")
cat("  Local coupling (adjacent wards): ", USE_LOCAL_COUPLING, "\n")
cat("  Rebalancing method: ", REBALANCE_METHOD, "\n")


# Check if local coupling is available
has_coupling <- exists("W_adj") && !is.null(W_adj)

# ==============================================================================
# Override Aedes recruitment + VPH from posterior median (standing-water model)
# ==============================================================================
# Replaces the legacy pre-scaled MAT_A_BURNIN / MAT_A_EXP (built from rain_14d)
# with the standing-water version: daily rainfall -> exponential water stock
# with per-district half-life -> ramp recruitment scaling. Posterior median
# parameters are pulled from output/abc_posterior.rds. Culex inputs are left
# unchanged (rain_14d driver).
#
# Falls back to the legacy rain_14d matrices if no posterior is available --
# useful for first-time exploration, but the calibrated model uses the
# standing-water branch.

`%||%` <- function(x, y) if (is.null(x)) y else x
POSTERIOR_PATH <- Sys.getenv("ABC_POSTERIOR_PATH",
                             unset = "output/abc_posterior.rds")
# Optional: pick a specific sample from the accepted_grid by its sample_id.
# Set SAMPLE_ID to use that sample's parameters instead of the posterior
# median -- useful for sanity-checking a known "best" sample, since the median
# of a ridge-shaped posterior can fall in a worse-fitting region.
SAMPLE_ID <- Sys.getenv("SAMPLE_ID", unset = "")
SAMPLE_ID <- if (nzchar(SAMPLE_ID)) as.integer(SAMPLE_ID) else NULL

if (file.exists(POSTERIOR_PATH)) {
  abc <- readRDS(POSTERIOR_PATH)
  fixed <- abc$fixed_params %||% list()

  needed <- c("vph_a_kili", "vph_c_kili", "vph_a_babati", "vph_c_babati",
              "half_life_kili", "half_life_babati")

  if (!is.null(SAMPLE_ID)) {
    # Pick a specific sample by sample_id from the accepted_grid
    grid <- abc$accepted_grid %||% stop(
      "accepted_grid not present in posterior file; cannot resolve SAMPLE_ID."
    )
    row <- grid[grid$sample_id == SAMPLE_ID, , drop = FALSE]
    if (nrow(row) != 1) {
      stop(sprintf("SAMPLE_ID %d not found in accepted_grid (got %d rows).",
                   SAMPLE_ID, nrow(row)))
    }
    missing <- setdiff(needed, names(row))
    if (length(missing) > 0) {
      stop("accepted_grid is missing required parameters: ",
           paste(missing, collapse = ", "))
    }
    pmed <- unlist(row[, needed, drop = FALSE])
    cat("\nUsing SAMPLE_ID =", SAMPLE_ID, "from", POSTERIOR_PATH, "\n")
  } else {
    pmed <- apply(as.data.frame(abc$posterior), 2, median)
    missing <- setdiff(needed, names(pmed))
    if (length(missing) > 0) {
      stop("Posterior is missing required parameters: ",
           paste(missing, collapse = ", "),
           ".\nRebuild output/abc_posterior.rds from the standing-water calibration.")
    }
    cat("\nLoaded posterior median from", POSTERIOR_PATH, "\n")
  }
  print(round(pmed[needed], 2))

  # Per-ward half-lives, expanded by district membership
  ward_half_life <- ifelse(
    ward_district[ward_ids] == "babati",
    pmed[["half_life_babati"]], pmed[["half_life_kili"]]
  )
  names(ward_half_life) <- ward_ids

  W_MAX_FIXED <- fixed$w_max  %||% 200
  BASE_A      <- fixed$base_frac_A %||% 0.05

  aedes_mats <- build_aedes_water_driver(
    rain_daily_burn     = MAT_RAIN_DAILY_BURNIN,
    rain_daily_exp      = MAT_RAIN_DAILY_EXP,
    half_lives_per_ward = ward_half_life,
    w_max               = W_MAX_FIXED,
    base_frac_A         = BASE_A
  )
  MAT_A_BURNIN <- aedes_mats$MAT_A_BURN
  MAT_A_EXP    <- aedes_mats$MAT_A_EXP

  # Per-ward VPH from posterior median; override params_static carrying capacities
  vph_a <- ifelse(ward_district[ward_ids] == "babati",
                  pmed[["vph_a_babati"]], pmed[["vph_a_kili"]])
  vph_c <- ifelse(ward_district[ward_ids] == "babati",
                  pmed[["vph_c_babati"]], pmed[["vph_c_kili"]])
  names(vph_a) <- names(vph_c) <- ward_ids

  params_static$lam_A_max <- params_static$mu_A * vph_a
  params_static$lam_C_max <- params_static$mu_C * vph_c
  params_static$eps_local <- fixed$eps_local %||% params_static$eps_local
  params_static$omega_h   <- fixed$omega_h   %||% params_static$omega_h

  cat("  VPH (per-ward, by district):",
      "Kili A=", pmed[["vph_a_kili"]],   "C=", pmed[["vph_c_kili"]],   "|",
      "Babati A=", pmed[["vph_a_babati"]], "C=", pmed[["vph_c_babati"]], "\n")
  cat("  Aedes half-lives: kili =", pmed[["half_life_kili"]],
      "d, babati =", pmed[["half_life_babati"]], "d\n")
  cat("  Saturation cap W_max =", W_MAX_FIXED, "\n")
} else {
  message("\nNo posterior at ", POSTERIOR_PATH,
          " -- falling back to legacy rain_14d-based MAT_A_* matrices.\n",
          "Run run_abc_posterior.R first for the calibrated standing-water model.")
}

# ==============================================================================
# Climate overview (for hindcast validation)
# ==============================================================================
# shows whether climate data captures known anomalies (e.g. 2006 El Niño rainfall).

cat("\n", strrep("=", 70), "\n")
cat("Climate overview\n")
cat(strrep("=", 70), "\n")

p_climate_overview <- plot_climate_overview(
  ward_climate = ward_climate,
  ward_district = ward_district,
  ndvi_df = if (exists("ward_ndvi") && !is.null(ward_ndvi)) ward_ndvi else NULL,
  outbreak_years = c(2006, 2007)
)
print(p_climate_overview)
# ggsave("figures/climate_overview_hindcast.png", p_climate_overview,
#        width = 10, height = 10, dpi = 300)

# ==============================================================================
# Visualise climate response curves
# ==============================================================================

sample_ward <- ward_ids[1]
message("  Sample ward for plots: ", sample_ward)

rainfall_plots <- plot_rainfall_curves(
  climate_df = ward_climate,
  sample_ward = sample_ward
)
rainfall_plots$combined
# ggsave("figures/rainfall_dependent_recruitment.png", width = 8, height = 8, dpi = 300)

plot_thermal_curves()
# ggsave("figures/thermal_performance_curves.png", width = 8, height = 6, dpi = 300)

# NOTE: habitat is currently uniform for all wards (no spatial multipliers).
# This script therefore runs a single baseline simulation. Counterfactual
# comparisons are handled in dedicated scripts such as run_figures.R.


# ==============================================================================
# Initial state
# ==============================================================================

set.seed(422)
y0 <- make_initial_state(ward_ids, N_star, params_static)
message("\nInitial state: ", length(y0), " variables (", n_wards, " wards x 13 compartments)")

# ==============================================================================
# Helper: run burn-in + experiment
# ==============================================================================

# Run the baseline model: burn-in to equilibrium, then the 2023-2024
# experiment period. Burn-in uses cycled climate so the model approaches a
# periodic endemic state before seeing the experiment climate.
run_scenario <- function(scenario_name,
                         MAT_HAB_A_BURNIN, MAT_HAB_C_BURNIN,
                         MAT_HAB_A_EXP, MAT_HAB_C_EXP) {
  cat("\n  Burn-in (", EXP_BURNIN_YEARS, " years) -> endemic state\n")

  burnin_dates <- seq.Date(as.Date("2000-01-01"), by = "day", length.out = T_BURNIN)

  out_burnin <- run_ode(
    MAT_A = MAT_A_BURNIN, MAT_C = MAT_C_BURNIN,
    MAT_TEMP = MAT_TEMP_BURNIN,
    dates = burnin_dates,
    y_init = y0,
    params = params_static,
    use_movement = USE_MOVEMENT,
    move_mats = move_mats,
    move_freq = MOVE_FREQ,
    move_step_days = MOVE_STEP_DAYS,
    move_scale_to_time = MOVE_SCALE_TO_TIME,
    move_mult = MOVE_MULT,
    ward_ids = ward_ids,
    N_star = N_star,
    rebalance_hosts = REBALANCE_HOSTS,
    rebalance_method = REBALANCE_METHOD,
    use_temp_dependence = USE_TEMP_DEPENDENCE,
    use_local_coupling = has_coupling && USE_LOCAL_COUPLING,
    W_adj = if (has_coupling) W_adj else NULL,
    adj_deg = if (has_coupling) adj_deg else NULL,
    MAT_HABITAT_A = MAT_HAB_A_BURNIN,
    MAT_HABITAT_C = MAT_HAB_C_BURNIN
  )

  state_cols <- setdiff(names(out_burnin), c("date", "time"))
  endemic_state <- as.numeric(out_burnin[nrow(out_burnin), state_cols])
  names(endemic_state) <- state_cols
  endemic_state <- reset_cumulative(endemic_state, ward_ids)

  endemic_prev <- calc_ward_endemic_prevalence(endemic_state, ward_ids)
  cat("  Mean ward prevalence:", sprintf("%.3f%%", mean(endemic_prev$prevalence) * 100), "\n")
  cat("  Mean seroprevalence:", sprintf("%.1f%%", mean(endemic_prev$seroprevalence) * 100), "\n")

  cat(
    "\n  Running experiment (", as.character(EXP_START), "to",
    as.character(EXP_END), ")...\n"
  )

  out_exp <- run_ode(
    MAT_A = MAT_A_EXP, MAT_C = MAT_C_EXP,
    MAT_TEMP = MAT_TEMP_EXP,
    dates = exp_dates,
    y_init = endemic_state,
    params = params_static,
    use_movement = USE_MOVEMENT,
    move_mats = move_mats,
    move_freq = MOVE_FREQ,
    move_step_days = MOVE_STEP_DAYS,
    move_scale_to_time = MOVE_SCALE_TO_TIME,
    move_mult = MOVE_MULT,
    ward_ids = ward_ids,
    N_star = N_star,
    rebalance_hosts = REBALANCE_HOSTS,
    rebalance_method = REBALANCE_METHOD,
    use_temp_dependence = USE_TEMP_DEPENDENCE,
    use_local_coupling = has_coupling && USE_LOCAL_COUPLING,
    W_adj = if (has_coupling) W_adj else NULL,
    adj_deg = if (has_coupling) adj_deg else NULL,
    MAT_HABITAT_A = MAT_HAB_A_EXP,
    MAT_HABITAT_C = MAT_HAB_C_EXP
  )

  metrics <- calc_global_metrics_from_out(out_exp, ward_ids)
  cat("  Peak prevalence:", sprintf("%.4f%%", metrics$peak_prev * 100), "\n")
  cat("  Cumulative incidence:", format(round(metrics$cum_incidence), big.mark = ","), "\n")

  list(
    out_burnin = out_burnin,
    out_exp = out_exp,
    endemic_state = endemic_state,
    endemic_prev = endemic_prev,
    metrics = metrics,
    scenario = scenario_name
  )
}

# ==============================================================================
# Baseline simulation
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Baseline simulation\n")
cat(strrep("=", 70), "\n")

res_baseline <- run_scenario(
  "Baseline",
  MAT_HAB_A_BURNIN = MAT_HABITAT_A_BURNIN,
  MAT_HAB_C_BURNIN = MAT_HABITAT_C_BURNIN,
  MAT_HAB_A_EXP = MAT_HABITAT_A_EXP,
  MAT_HAB_C_EXP = MAT_HABITAT_C_EXP
)

# ==============================================================================
# Infection-probability calculation
# ==============================================================================
# Window-specific infection probability = 1 - exp(-cumulative hazard).
# Computed from the per-genus cumulative hazard trackers in the ODE output.
# We compare the 12-month pre-baseline 2023 window with the 8-month 2024
# post-baseline window used in the v3 calibration.

cat("\n", strrep("=", 70), "\n")
cat("Force of infection\n")
cat(strrep("=", 70), "\n")

# v3 survey-aligned windows.
d23_start <- as.Date("2022-09-14")
d23_end   <- as.Date("2023-09-14")
d24_start <- as.Date("2023-09-27")
d24_end   <- as.Date("2024-05-27")

foi_all <- bind_rows(
  summarise_annual_foi_from_Hh(res_baseline$out_exp, ward_ids, d23_start, d23_end) |>
    mutate(year = "2023"),
  summarise_annual_foi_from_Hh(res_baseline$out_exp, ward_ids, d24_start, d24_end) |>
    mutate(year = "2024")
) |>
  mutate(district = ward_district[ward])

ward_coords <- get_ward_coords(ward_info)
if (!is.null(ward_coords)) {
  foi_all <- foi_all |> left_join(ward_coords, by = "ward")
}

# ==============================================================================
# Posterior-median validation: model vs observed at ward and district level
# ==============================================================================
# Produces output/fig_validation_median.pdf -- a quick sanity check that the
# posterior median reproduces the 2023 and 2024 spatial FOI patterns. The
# full posterior-uncertainty version (with credible intervals) is built by
# run_figures.R as fig0_validation.pdf.

cat("\n", strrep("=", 70), "\n")
cat("Posterior-median validation (model vs observed)\n")
cat(strrep("=", 70), "\n")

cal_targets <- build_calibration_targets(N_star = N_star)
surveyed_wards <- cal_targets$surveyed_wards
ws <- cal_targets$ws

# Ward-level: model probability vs observed probability target, by year.
# 2023: model `1-exp(-H over 12mo pre-baseline)` ↔ survey `seroconv_hat_2023`.
# 2024: model `1-exp(-H over 8mo)` ↔ survey `seroconv_hat_2024_8mo`.
# All on probability scale (v3).
ward_val <- foi_all |>
  filter(ward %in% surveyed_wards) |>
  select(ward, district, year, model_foi = annual_foi_prob) |>
  left_join(
    ws |> dplyr::select(ward = clean_name,
                        seroconv_hat_2023, seroconv_hat_2024_8mo),
    by = "ward"
  ) |>
  mutate(obs_foi = ifelse(year == "2023",
                          seroconv_hat_2023, seroconv_hat_2024_8mo)) |>
  filter(!is.na(obs_foi), obs_foi > 0)

# Pearson correlation per year (printed in the scatter panels)
corr_tab <- ward_val |>
  group_by(year) |>
  summarise(corr = cor(model_foi, obs_foi, use = "complete.obs"),
            n = n(), .groups = "drop")
cat("\nWard-level correlation (model vs observed):\n")
print(corr_tab)

# Panel A/B: ward scatter, model vs observed (one facet per year)
p_scatter <- ward_val |>
  ggplot(aes(x = obs_foi * 100, y = model_foi * 100, colour = district)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  facet_wrap(~year, scales = "free") +
  geom_text(
    data = corr_tab,
    aes(x = Inf, y = -Inf, label = sprintf("r = %.2f", corr)),
    colour = "black", hjust = 1.1, vjust = -0.5,
    size = 3.6, inherit.aes = FALSE
  ) +
  theme_bw(base_size = 11) +
  labs(
    x = "Observed infection probability (%)",
    y = "Posterior-median model infection probability (%)",
    colour = "District",
    title = "Ward-level infection probability: posterior median vs observed"
  )

# Panel C: district mean comparison (simple mean across surveyed wards,
# matching the cluster calibration estimator).
district_compare <- ward_val |>
  group_by(year, district) |>
  summarise(
    model_mean = mean(model_foi),
    obs_mean   = mean(obs_foi),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    c(model_mean, obs_mean),
    names_to = "source", values_to = "foi"
  ) |>
  mutate(source = ifelse(source == "model_mean",
                         "Model (posterior median)", "Observed"))

cat("\nDistrict mean infection probability (%), model vs observed (simple ward-mean):\n")
print(district_compare |>
        tidyr::pivot_wider(names_from = source, values_from = foi) |>
        mutate(across(where(is.numeric), ~ round(.x * 100, 2))))

p_district <- district_compare |>
  ggplot(aes(x = district, y = foi * 100, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~year, scales = "free_y") +
  scale_fill_manual(values = c(
    "Model (posterior median)" = "steelblue",
    "Observed"                 = "grey30"
  )) +
  theme_bw(base_size = 11) +
  labs(
    x = NULL, y = "Infection probability (%)", fill = NULL,
    title = "District mean infection probability: model vs observed (surveyed wards, simple mean)"
  )

# Panel D: ward-level cumulative seroprevalence distribution, 2023 only.
# 2024 dropped here because Paul's 2024 estimand is 8-month seroconversion
# (already shown in the FOI panel), not cumulative seroprev. Comparing
# model cumulative seroprev at the 2024 survey midpoint to Paul's
# seroconv_2024_8mo would mix scales.
sero_day_23 <- which(exp_dates == as.Date("2023-09-14"))

model_sero_long <- extract_ward_seroprevalence(
    res_baseline$out_exp, ward_ids, sero_day_23
  ) |>
  filter(ward_id %in% surveyed_wards) |>
  mutate(district = ward_district[ward_id], year = "2023") |>
  transmute(district, year, sero = seroprev,
            source = "Model (posterior median)")

obs_sero_long <- ws |>
  filter(clean_name %in% surveyed_wards) |>
  dplyr::select(district, seroprev_2023) |>
  mutate(year = "2023", source = "Observed", sero = seroprev_2023) |>
  dplyr::select(district, year, sero, source) |>
  filter(!is.na(sero), sero > 0)

sero_long <- bind_rows(model_sero_long, obs_sero_long)

p_sero <- sero_long |>
  ggplot(aes(x = district, y = sero * 100, fill = source)) +
  geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA,
               position = position_dodge(0.7)) +
  geom_point(aes(colour = source),
             position = position_jitterdodge(
               dodge.width = 0.7, jitter.width = 0.10
             ),
             size = 1.3, alpha = 0.75) +
  facet_wrap(~year, scales = "free_y") +
  scale_fill_manual(values = c(
    "Model (posterior median)" = "steelblue",
    "Observed"                 = "grey30"
  )) +
  scale_colour_manual(values = c(
    "Model (posterior median)" = "steelblue",
    "Observed"                 = "grey30"
  ), guide = "none") +
  theme_bw(base_size = 11) +
  labs(
    x = NULL, y = "Seroprevalence (%)", fill = NULL,
    title = "Ward seroprevalence by district: model vs observed (distribution over wards)"
  )

# Quick console summary for sero
sero_summary <- sero_long |>
  group_by(year, district, source) |>
  summarise(median = median(sero, na.rm = TRUE),
            q25 = quantile(sero, 0.25, na.rm = TRUE),
            q75 = quantile(sero, 0.75, na.rm = TRUE),
            n = n(), .groups = "drop") |>
  mutate(across(c(median, q25, q75), ~ round(.x * 100, 2)))
cat("\nWard seroprevalence (%) summary, model vs observed:\n")
print(sero_summary)

p_validation_median <- (p_scatter / p_district / p_sero) +
  plot_layout(heights = c(1.2, 1, 1))
ggsave("output/fig_validation_median.pdf",
       p_validation_median, width = 11, height = 12)
cat("\nSaved validation plot -> output/fig_validation_median.pdf\n")

# ==============================================================================
# Diagnostic plots
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Generating diagnostic plots\n")
cat(strrep("=", 70), "\n")

# --- Vectors per host ---
cat("\n  Vectors per host...\n")

vph_ts <- extract_regional_ts(res_baseline$out_exp, "Baseline") |>
  select(date, Aedes, Culex) |>
  pivot_longer(cols = c(Aedes, Culex), names_to = "species", values_to = "vph")

p_vph <- vph_ts |>
  ggplot(aes(x = date, y = vph, color = species)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  scale_color_manual(values = cols_vec) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    title = "Vectors per host", x = NULL, y = "Vectors per host",
    color = "Species"
  )

print(p_vph)
# ggsave("figures/vectors_per_host.png", p_vph, width = 10, height = 5, dpi = 300)

# --- Spatial FOI maps ---
cat("\n  FOI spatial maps...\n")

foi_2023_range <- foi_all |>
  filter(year == "2023") |>
  pull(annual_foi_prob) |>
  range(na.rm = TRUE)
foi_2024_range <- foi_all |>
  filter(year == "2024") |>
  pull(annual_foi_prob) |>
  range(na.rm = TRUE)
foi_range <- range(foi_all$annual_foi_prob, na.rm = TRUE)

cat("  2023 infection-probability range:",
    sprintf("%.1f%% - %.1f%%", foi_2023_range[1] * 100, foi_2023_range[2] * 100),
    "\n")
cat("  2024 infection-probability range:",
    sprintf("%.1f%% - %.1f%%", foi_2024_range[1] * 100, foi_2024_range[2] * 100),
    "\n")

make_foi_map <- function(data, title_str, color_limits) {
  data |>
    ggplot(aes(x = lon, y = lat, color = annual_foi_prob, shape = district)) +
    geom_point(size = 3.5, alpha = 0.9) +
    scale_color_viridis_c(
      option = "H", limits = color_limits,
      labels = percent_format(), name = "Infection\nprobability"
    ) +
    coord_equal() +
    theme_bw() +
    theme(legend.position = "right", plot.title = element_text(size = 11)) +
    labs(title = title_str, x = NULL, y = NULL)
}

p_foi_maps <- (
  make_foi_map(
    foi_all |> filter(year == "2023"),
    "2023", foi_range
  ) |
    make_foi_map(
      foi_all |> filter(year == "2024"),
      "2024", foi_range
    )
) +
  plot_annotation(title = "Window-specific infection probability")

print(p_foi_maps)
# ggsave("figures/foi_spatial_maps.png", p_foi_maps, width = 12, height = 6, dpi = 300)

# ==============================================================================
# Summary statistics
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Summary statistics\n")
cat(strrep("=", 70), "\n")

foi_summary <- foi_all |>
  group_by(year) |>
  summarise(
    mean_foi = mean(annual_foi_prob, na.rm = TRUE),
    sd_foi = sd(annual_foi_prob, na.rm = TRUE),
    cv_foi = sd_foi / mean_foi,
    .groups = "drop"
  )

cat("\nWindow-specific infection-probability summary:\n")
cat("  Year | Mean prob | SD    | CV (heterogeneity)\n")
cat("  -----|-----------|-------|-------------------\n")
for (i in 1:nrow(foi_summary)) {
  cat(sprintf(
    "  %4s | %7.1f%% | %5.1f%% | %.2f\n",
    foi_summary$year[i],
    foi_summary$mean_foi[i] * 100, foi_summary$sd_foi[i] * 100,
    foi_summary$cv_foi[i]
  ))
}

cat("\nKey insight: CV measures spatial heterogeneity in infection risk.\n")

# ==============================================================================
# Standard plots
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Standard plots\n")
cat(strrep("=", 70), "\n")

# District time series (burn-in)
regional_ts_temp <- extract_regional_ts(res_baseline$out_burnin, "burn-in")
burnin_dates <- seq.Date(as.Date("2000-01-01"), by = "day", length.out = T_BURNIN)
burnin_end_date_long <- max(burnin_dates)

district_ts_temp <- res_baseline$out_burnin |>
  select(
    time, date, starts_with("Sh__"), starts_with("Eh__"),
    starts_with("Ih__"), starts_with("Rh__")
  ) |>
  pivot_longer(-c(time, date), names_to = c("compartment", "ward"), names_sep = "__") |>
  mutate(district = ward_district[ward]) |>
  group_by(time, date, district, compartment) |>
  summarise(total = sum(value), .groups = "drop") |>
  pivot_wider(names_from = compartment, values_from = total) |>
  mutate(Nh = Sh + Eh + Ih + Rh, seroprevalence = Rh / pmax(Nh, 1))

p1 <- district_ts_temp |>
  filter(!is.na(district)) |>
  ggplot(aes(x = date, y = seroprevalence, color = district)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = burnin_end_date_long, linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = cols_district,
    labels = c("babati" = "Babati", "hai" = "Hai", "moshi" = "Moshi")
  ) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  labs(title = "Seroprevalence (burn-in)", x = NULL, y = "Seroprevalence", color = "District")

p2 <- regional_ts_temp |>
  filter(year(date) >= max(year(date))) |>
  select(date, Aedes, Culex) |>
  pivot_longer(-date, names_to = "species", values_to = "count") |>
  ggplot(aes(x = date, y = count, color = species)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = cols_vec) +
  theme_bw() +
  labs(
    title = "Vectors per host (last year of burn-in)", x = NULL,
    y = "Vectors/host", color = NULL
  )

p3 <- regional_ts_temp |>
  ggplot(aes(x = date, y = prevalence)) +
  geom_line(color = "#e74c3c", linewidth = 0.8) +
  geom_vline(xintercept = burnin_end_date_long, linetype = "dashed", color = "gray40") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  labs(title = "Regional prevalence (burn-in)", x = NULL, y = "Prevalence")

print((p1 / p3) / p2)
# ggsave("figures/endemic_burnin_timeseries.png", width = 8, height = 8, dpi = 300)

# Baseline seroprevalence (2023-2024)
district_ts <- build_district_ts(res_baseline$out_exp, exp_dates, ward_ids, ward_district)

p_baseline <- district_ts |>
  ggplot(aes(x = date, y = sero, color = district)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = cols_district,
    labels = c("babati" = "Babati", "hai" = "Hai", "moshi" = "Moshi")
  ) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  labs(
    title = "Baseline seroprevalence (2023-2024)", x = NULL,
    y = "Seroprevalence", color = "District"
  )

print(p_baseline)
# ggsave("figures/baseline_seroprevalence_2023_2024.png", width = 8, height = 6, dpi = 300)

# FOI boxplots
foi_plot <- foi_all |>
  mutate(year_label = ifelse(year == "2023", "2023 (Baseline)", "2024 (El Nino)"))

p_box_prob <- foi_plot |>
  ggplot(aes(x = district, y = annual_foi_prob, fill = district)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2),
    size = 1.5, alpha = 0.4, aes(color = district)
  ) +
  facet_wrap(~year_label, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  labs(
    title = "Annual risk of infection",
    x = NULL, y = "Annual probability (%)"
  )

print(p_box_prob)
# ggsave("figures/foi_annual_boxplots.png", width = 8, height = 6, dpi = 300)

# FOI maps
p_map_23 <- foi_plot |>
  filter(year == "2023") |>
  ggplot(aes(x = lon, y = lat, color = annual_foi_prob, shape = district)) +
  geom_point(size = 3.5, alpha = 1) +
  coord_equal() +
  theme_bw() +
  scale_color_viridis_c(option = "H", direction = 1, labels = percent_format()) +
  labs(title = "2023 (Baseline)", x = NULL, y = NULL, color = "FOI") +
  theme(legend.position = "right")

p_map_24 <- foi_plot |>
  filter(year == "2024") |>
  ggplot(aes(x = lon, y = lat, color = annual_foi_prob, shape = district)) +
  geom_point(size = 3.5, alpha = 1) +
  coord_equal() +
  theme_bw() +
  scale_color_viridis_c(option = "H", direction = 1, labels = percent_format()) +
  labs(title = "2024 (El Nino)", x = NULL, y = NULL, color = "FOI") +
  theme(legend.position = "right")

print((p_map_23 / p_map_24) +
  plot_annotation(title = "Spatial distribution of infection risk"))
# ggsave("figures/foi_annual_maps.png", width = 8, height = 8, dpi = 300)


# ==============================================================================
# HINDCAST: 2000-2024 with actual climate (optional, see run_hindcast.R)
# ==============================================================================
# Alternative to run_hindcast.R — uses a simpler burn-in (cycling 2000-2004
# climate only) rather than the full 3-cycle strategy in run_hindcast.R.
#
# Validates model against documented RVF outbreaks:
#   - 2006-2007: Major East Africa outbreak (primary target)
#   - Non-outbreak years should show low incidence
#
# Key finding: 2007 does NOT stand out as a climate anomaly, so the model
# won't produce a spike unless case importation is modelled explicitly.
#
# Approach:
#   Phase 1: Burn-in with cycled 2000-2004 climate -> endemic equilibrium
#   Phase 2: Forward run with actual sequential climate 2000-2024
#   Phase 3: Compare predicted incidence peaks to documented outbreaks

if (EXPERIMENT_MODE == "hindcast" &&
  exists("MAT_A_HINDCAST") && !is.null(MAT_A_HINDCAST)) {
  cat("\n", strrep("=", 70), "\n")
  cat("Hindcast: 2000-2024 forward run\n")
  cat(strrep("=", 70), "\n")

  # --------------------------------------------------------------------------
  # Phase 1: Burn-in to endemic equilibrium
  # --------------------------------------------------------------------------
  # Uses 2000-2004 climate cycled for HINDCAST_BURNIN_YEARS.
  # The model never sees the 2006 rainfall anomaly during burn-in.

  cat(
    "\n  Phase 1: Burn-in (", HINDCAST_BURNIN_YEARS,
    " years, cycling ", paste(range(HINDCAST_BURNIN_CLIMATE_YEARS), collapse = "-"),
    " climate)\n"
  )

  hindcast_burnin_dates <- seq.Date(
    as.Date("1990-01-01"),
    by = "day", length.out = T_HINDCAST_BURNIN
  )

  out_hindcast_burnin <- run_ode(
    MAT_A = MAT_A_HINDCAST_BURNIN,
    MAT_C = MAT_C_HINDCAST_BURNIN,
    MAT_TEMP = MAT_TEMP_HINDCAST_BURNIN,
    dates = hindcast_burnin_dates,
    y_init = y0,
    params = params_static,
    use_movement = USE_MOVEMENT,
    move_mats = move_mats,
    move_freq = MOVE_FREQ,
    move_step_days = MOVE_STEP_DAYS,
    move_scale_to_time = MOVE_SCALE_TO_TIME,
    move_mult = MOVE_MULT,
    ward_ids = ward_ids,
    N_star = N_star,
    rebalance_hosts = REBALANCE_HOSTS,
    rebalance_method = REBALANCE_METHOD,
    use_temp_dependence = USE_TEMP_DEPENDENCE,
    use_local_coupling = has_coupling && USE_LOCAL_COUPLING,
    W_adj = if (has_coupling) W_adj else NULL,
    adj_deg = if (has_coupling) adj_deg else NULL,
    MAT_HABITAT_A = MAT_HABITAT_A_HINDCAST,
    MAT_HABITAT_C = MAT_HABITAT_C_HINDCAST
  )

  # Extract endemic state
  state_cols_hc <- setdiff(names(out_hindcast_burnin), c("date", "time"))
  ENDEMIC_STATE_HINDCAST <- as.numeric(
    out_hindcast_burnin[nrow(out_hindcast_burnin), state_cols_hc]
  )
  names(ENDEMIC_STATE_HINDCAST) <- state_cols_hc
  ENDEMIC_STATE_HINDCAST <- reset_cumulative(ENDEMIC_STATE_HINDCAST, ward_ids)

  endemic_prev_hc <- calc_ward_endemic_prevalence(ENDEMIC_STATE_HINDCAST, ward_ids)
  cat("  Endemic prevalence:", sprintf("%.3f%%", mean(endemic_prev_hc$prevalence) * 100), "\n")
  cat("  Endemic seroprevalence:", sprintf("%.1f%%", mean(endemic_prev_hc$seroprevalence) * 100), "\n")

  # --------------------------------------------------------------------------
  # Phase 2: Forward run with actual climate (2000-2024)
  # --------------------------------------------------------------------------
  # This is the key validation step. The model encounters 2006 anomalous
  # rainfall for the first time. We check whether incidence spikes.

  cat(
    "\n  Phase 2: Forward run (", as.character(HINDCAST_START), "to",
    as.character(HINDCAST_END), ",", length(hindcast_dates), "days)\n"
  )

  out_hindcast <- run_ode(
    MAT_A = MAT_A_HINDCAST,
    MAT_C = MAT_C_HINDCAST,
    MAT_TEMP = MAT_TEMP_HINDCAST,
    dates = hindcast_dates,
    y_init = ENDEMIC_STATE_HINDCAST,
    params = params_static,
    use_movement = USE_MOVEMENT,
    move_mats = move_mats,
    move_freq = MOVE_FREQ,
    move_step_days = MOVE_STEP_DAYS,
    move_scale_to_time = MOVE_SCALE_TO_TIME,
    move_mult = MOVE_MULT,
    ward_ids = ward_ids,
    N_star = N_star,
    rebalance_hosts = REBALANCE_HOSTS,
    rebalance_method = REBALANCE_METHOD,
    use_temp_dependence = USE_TEMP_DEPENDENCE,
    use_local_coupling = has_coupling && USE_LOCAL_COUPLING,
    W_adj = if (has_coupling) W_adj else NULL,
    adj_deg = if (has_coupling) adj_deg else NULL,
    MAT_HABITAT_A = MAT_HABITAT_A_HINDCAST,
    MAT_HABITAT_C = MAT_HABITAT_C_HINDCAST
  )

  cat("  Forward run complete.\n")

  # --------------------------------------------------------------------------
  # Extract time series for plotting
  # --------------------------------------------------------------------------

  # Regional aggregation
  hindcast_regional <- extract_regional_ts(out_hindcast, "hindcast")

  # District-level daily incidence from cumulative tracker
  hindcast_Ch <- extract_comp_mat(out_hindcast, "Ch", ward_ids)
  hindcast_daily_inc <- apply(hindcast_Ch, 2, function(x) pmax(0, c(x[1], diff(x))))

  # Aggregate daily incidence to district level
  hindcast_inc_district <- tibble(
    date = hindcast_dates,
    time = seq_along(hindcast_dates)
  ) |>
    bind_cols(as_tibble(hindcast_daily_inc)) |>
    pivot_longer(-c(date, time), names_to = "ward", values_to = "daily_inc") |>
    mutate(
      ward = str_remove(ward, "^Ch__"),
      district = ward_district[ward]
    ) |>
    filter(!is.na(district)) |>
    group_by(date, district) |>
    summarise(daily_inc = sum(daily_inc, na.rm = TRUE), .groups = "drop")

  # Monthly incidence (smoother for visualization)
  hindcast_monthly_inc <- hindcast_inc_district |>
    mutate(year_month = floor_date(date, "month")) |>
    group_by(year_month, district) |>
    summarise(monthly_inc = sum(daily_inc, na.rm = TRUE), .groups = "drop")

  # Annual incidence per district
  hindcast_annual_inc <- hindcast_inc_district |>
    mutate(year = year(date)) |>
    group_by(year, district) |>
    summarise(annual_inc = sum(daily_inc, na.rm = TRUE), .groups = "drop")

  # Annual FOI per ward per year
  hindcast_years <- sort(unique(year(hindcast_dates)))
  hindcast_foi_list <- list()

  for (yr in hindcast_years) {
    yr_start <- as.Date(sprintf("%d-01-01", yr))
    yr_end <- as.Date(sprintf("%d-12-31", yr))

    # Check dates exist in output
    if (yr_start < min(hindcast_dates) || yr_end > max(hindcast_dates)) next

    foi_yr <- tryCatch(
      summarise_annual_foi_from_Hh(out_hindcast, ward_ids, yr_start, yr_end) |>
        mutate(year = yr),
      error = function(e) NULL
    )
    if (!is.null(foi_yr)) hindcast_foi_list[[length(hindcast_foi_list) + 1]] <- foi_yr
  }

  hindcast_foi_annual <- bind_rows(hindcast_foi_list) |>
    mutate(district = ward_district[ward])

  # --------------------------------------------------------------------------
  # PLOT H1: Monthly incidence by district (full time series)
  # --------------------------------------------------------------------------
  cat("\n  Generating hindcast plots...\n")

  # Known outbreak periods for reference
  outbreak_periods <- tibble(
    xmin = as.Date(c("2006-11-01", "2014-01-01", "2018-01-01")),
    xmax = as.Date(c("2007-06-01", "2015-06-01", "2018-06-01")),
    label = c("2006-07 outbreak", "2014-15", "2018")
  )

  p_h1_monthly <- hindcast_monthly_inc |>
    ggplot(aes(x = year_month, y = monthly_inc, color = district)) +
    geom_rect(
      data = outbreak_periods,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "red", alpha = 0.08
    ) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(
      values = cols_district,
      labels = c("babati" = "Babati", "hai" = "Hai", "moshi" = "Moshi")
    ) +
    theme_bw() +
    labs(
      title = "Hindcast: Monthly incidence by district",
      subtitle = "Red shading = documented RVF outbreak periods",
      x = NULL, y = "Monthly new infections", color = "District"
    )

  # --------------------------------------------------------------------------
  # PLOT H2: Annual incidence bar chart
  # --------------------------------------------------------------------------

  p_h2_annual <- hindcast_annual_inc |>
    ggplot(aes(x = year, y = annual_inc, fill = district)) +
    geom_col(position = "stack", alpha = 0.85) +
    geom_vline(
      xintercept = c(2006.5, 2007.5), linetype = "dashed",
      color = "red", alpha = 0.7
    ) +
    scale_fill_manual(
      values = cols_district,
      labels = c("babati" = "Babati", "hai" = "Hai", "moshi" = "Moshi")
    ) +
    theme_bw() +
    labs(
      title = "Hindcast: Annual incidence",
      subtitle = "Red dashed lines bracket 2006-2007 outbreak",
      x = "Year", y = "Total annual infections", fill = "District"
    )

  # --------------------------------------------------------------------------
  # PLOT H3: Regional prevalence and vectors per host
  # --------------------------------------------------------------------------

  p_h3_prev <- hindcast_regional |>
    ggplot(aes(x = date, y = prevalence)) +
    geom_rect(
      data = outbreak_periods,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "red", alpha = 0.08
    ) +
    geom_line(color = "#e74c3c", linewidth = 0.6) +
    scale_y_continuous(labels = percent_format()) +
    theme_bw() +
    labs(title = "Regional prevalence", x = NULL, y = "Prevalence")

  p_h3_vph <- hindcast_regional |>
    select(date, Aedes, Culex) |>
    pivot_longer(-date, names_to = "species", values_to = "vph") |>
    ggplot(aes(x = date, y = vph, color = species)) +
    geom_rect(
      data = outbreak_periods,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "red", alpha = 0.08
    ) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = cols_vec) +
    theme_bw() +
    labs(title = "Vectors per host", x = NULL, y = "VPH", color = NULL)

  # --------------------------------------------------------------------------
  # PLOT H4: Annual FOI heatmap (wards x years)
  # --------------------------------------------------------------------------

  p_h4_foi_heatmap <- hindcast_foi_annual |>
    ggplot(aes(x = year, y = reorder(ward, annual_foi_prob), fill = annual_foi_prob)) +
    geom_tile() +
    scale_fill_viridis_c(option = "H", labels = percent_format(), name = "Annual FOI") +
    geom_vline(xintercept = c(2006.5, 2007.5), color = "red", linewidth = 0.5) +
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(
      title = "Annual FOI by ward",
      subtitle = "Wards ordered by mean FOI; red lines = 2006-2007",
      x = "Year", y = paste0("Wards (n=", n_wards, ")")
    )

  # --------------------------------------------------------------------------
  # PLOT H5: Mean annual FOI across wards (are outbreak years elevated?)
  # --------------------------------------------------------------------------

  foi_annual_summary <- hindcast_foi_annual |>
    group_by(year) |>
    summarise(
      mean_foi = mean(annual_foi_prob, na.rm = TRUE),
      lo = quantile(annual_foi_prob, 0.25, na.rm = TRUE),
      hi = quantile(annual_foi_prob, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(outbreak = year %in% c(2006, 2007))

  p_h5_foi_annual <- foi_annual_summary |>
    ggplot(aes(x = year, y = mean_foi)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
    geom_line(color = "steelblue", linewidth = 0.9) +
    geom_point(aes(color = outbreak), size = 2.5) +
    scale_color_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "red"),
      labels = c("Non-outbreak", "Outbreak"),
      name = NULL
    ) +
    scale_y_continuous(labels = percent_format()) +
    theme_bw() +
    labs(
      title = "Mean annual FOI across wards",
      subtitle = "Ribbon = IQR across wards; red = documented outbreak years",
      x = "Year", y = "Annual infection probability"
    )

  # --------------------------------------------------------------------------
  # Combine and print
  # --------------------------------------------------------------------------

  # Main hindcast figure
  p_hindcast_main <- (p_h1_monthly / p_h2_annual) +
    plot_annotation(title = "Hindcast validation: 2000-2024")
  print(p_hindcast_main)
  # ggsave("figures/hindcast_incidence.png", p_hindcast_main, width = 12, height = 8, dpi = 300)

  # Dynamics figure
  p_hindcast_dynamics <- (p_h3_prev / p_h3_vph) +
    plot_annotation(title = "Hindcast: Regional dynamics")
  print(p_hindcast_dynamics)
  # ggsave("figures/hindcast_dynamics.png", p_hindcast_dynamics, width = 12, height = 6, dpi = 300)

  # FOI figure
  p_hindcast_foi <- (p_h5_foi_annual / p_h4_foi_heatmap) +
    plot_annotation(title = "Hindcast: Annual force of infection")
  print(p_hindcast_foi)
  # ggsave("figures/hindcast_foi.png", p_hindcast_foi, width = 12, height = 10, dpi = 300)

  # --------------------------------------------------------------------------
  # Summary table
  # --------------------------------------------------------------------------
  cat("\n  Hindcast annual incidence summary:\n")
  cat("  Year | Total inc. | Mean FOI  | Outbreak?\n")
  cat("  -----|------------|-----------|----------\n")

  annual_totals <- hindcast_annual_inc |>
    group_by(year) |>
    summarise(total = sum(annual_inc), .groups = "drop") |>
    left_join(foi_annual_summary |> select(year, mean_foi, outbreak), by = "year")

  for (i in 1:nrow(annual_totals)) {
    flag <- if (annual_totals$outbreak[i]) " ***" else ""
    cat(sprintf(
      "  %d | %10s | %8.1f%% | %s%s\n",
      annual_totals$year[i],
      format(round(annual_totals$total[i]), big.mark = ","),
      annual_totals$mean_foi[i] * 100,
      ifelse(annual_totals$outbreak[i], "YES", "no"),
      flag
    ))
  }

  # Outbreak vs non-outbreak discrimination
  foi_outbreak <- foi_annual_summary |>
    filter(outbreak) |>
    pull(mean_foi)
  foi_non_outbreak <- foi_annual_summary |>
    filter(!outbreak) |>
    pull(mean_foi)

  cat("\n  Outbreak years mean FOI:     ", sprintf("%.2f%%", mean(foi_outbreak) * 100), "\n")
  cat("  Non-outbreak years mean FOI: ", sprintf("%.2f%%", mean(foi_non_outbreak) * 100), "\n")
  cat(
    "  Fold increase:               ",
    sprintf("%.1fx", mean(foi_outbreak) / mean(foi_non_outbreak)), "\n"
  )
} else if (EXPERIMENT_MODE == "hindcast") {
  warning(
    "Hindcast mode selected but hindcast climate matrices not available.\n",
    "Check that climate data covers ", HINDCAST_START, " to ", HINDCAST_END, "."
  )
}
