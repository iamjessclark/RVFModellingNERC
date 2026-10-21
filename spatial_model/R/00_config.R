# ==============================================================================
# 00_config.R
#
# Configuration, parameters, and climate response functions for the
# ward-level SEIRS+SEI Rift Valley fever model (66 wards, 3 districts).
#
# Components:
#   - Rainfall-dependent recruitment (Aedes hill-shaped, Culex hump-shaped)
#   - Temperature-dependent vector biology (Fischer et al. 2013)
#   - Transmission competence (Aedes >> Culex per bite)
#   - District-level VPH carrying capacity (4 parameters)
#   - External introduction rate (background cattle trade)
# ==============================================================================

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
PATH_WARD_SUM <- "./ward_info/ward_summary_df.rds"
PATH_MOVEMENT <- "movement_data/sim_edge_list.rds"
PATH_CLIMATE <- "./climate_data/daily_climate_df.rds"
PATH_CONTIGUOUS <- "../SpatialData/contiguous.wards.csv"
PATH_NDVI <- "./climate_data/raw_ndvi_df.rds"

# ------------------------------------------------------------------------------
# Simulation switches
# ------------------------------------------------------------------------------
USE_MOVEMENT <- TRUE
USE_LOCAL_COUPLING <- TRUE
SIM_ID <- 1
REBALANCE_HOSTS <- TRUE
REBALANCE_METHOD <- "scale" # "scale" preserves immunity ratios; "susceptible" adds all to S

# Movement frequency and intensity controls
MOVE_FREQ <- "monthly" # or "fixed_step"
MOVE_STEP_DAYS <- 14 # 7 for weekly, 14 for biweekly
MOVE_SCALE_TO_TIME <- TRUE # keep annual movement constant when changing freq
MOVE_MULT <- 1 # sensitivity to movement intensity

# Temperature-dependent vector parameters
USE_TEMP_DEPENDENCE <- TRUE

# ------------------------------------------------------------------------------
# Hindcast configuration
# ------------------------------------------------------------------------------
# Burn-in strategy:
#   Cycle the full 2000-2024 climate data N_BURNIN_CYCLES times (2 times, ~75 years).
#   Forward run starts at 2000 from the second cycle.
#
# Why:
#   - The system sees realistic inter-annual variability during burn-in
#   - No arbitrary selection of "pre-outbreak" years
#   - 2006-2007 was not climatically exceptional, so no anomaly to hide from the burn-in

HINDCAST_START <- as.Date("2000-01-01")
HINDCAST_END <- as.Date("2024-12-31")
HINDCAST_T_DAYS <- as.integer(HINDCAST_END - HINDCAST_START) + 1

N_BURNIN_CYCLES <- 3 # Number of 2000-2024 cycles for burn-in

# --- Contemporary experiment (2022-09 to 2024-12) ---
# Start pushed back from 2023-01-01 (v2) to 2022-09-01 (v3) so the experiment
# window covers the survey-aligned calibration targets:
#   d23_start = 2022-09-14 (12 months pre-baseline survey)
#   d24_end   = 2024-05-27 (8 months ending at follow-up survey midpoint)
# The 2-week buffer before d23_start avoids edge effects when extracting
# cumulative hazard at the window boundary.
EXP_START <- as.Date("2022-09-01")
EXP_END <- as.Date("2024-12-31")
EXP_BURNIN_YEARS <- 40
EXP_T_DAYS <- as.integer(EXP_END - EXP_START) + 1

# For burn-in: which climate years to cycle for reaching endemic equilibrium?
# Uses a representative recent period that captures typical inter-annual
# variability without including the 2006-07 outbreak years.
# 2010-2022 = 13 years, cycled ~3 times over 40 years of burn-in.
BURNIN_CLIMATE_START <- as.Date("2010-01-01")
BURNIN_CLIMATE_END <- as.Date("2022-12-31")

# ------------------------------------------------------------------------------
# Plot colour palettes
# ------------------------------------------------------------------------------
cols_district <- c("babati" = "#e41a1c", "hai" = "#377eb8", "moshi" = "#4daf4a")
cols_host <- c(
  "S" = "#3498db",
  "E" = "#f39c12",
  "I" = "#e74c3c",
  "R" = "#27ae60"
)
cols_vec <- c("Aedes" = "firebrick", "Culex" = "forestgreen")

# ==============================================================================
# Rainfall-dependent recruitment functions
# ==============================================================================
# Aedes candidate curves are kept simple on purpose:
#   - baseline_hill: current working hill-shaped response
#   - early_flush_hill: same idea, but decline starts earlier
#   - briere_flush: smooth unimodal alternative to the piecewise hill
# These candidates let us test whether wetter Kilimanjaro wards need stronger
# suppression under excess rainfall, without hard-coding district-specific rules.
AEDES_CURVE_CANDIDATES <- c("baseline_hill", "early_flush_hill", "briere_flush")
AEDES_CURVE_CANDIDATE <- "baseline_hill"

get_aedes_curve_label <- function(candidate = AEDES_CURVE_CANDIDATE) {
  candidate <- match.arg(candidate, AEDES_CURVE_CANDIDATES)
  dplyr::case_match(
    candidate,
    "baseline_hill" ~ "Baseline hill",
    "early_flush_hill" ~ "Early-flush hill",
    "briere_flush" ~ "Briere-like flush"
  )
}

#' Aedes recruitment scaling (0 to 1): hill-shaped response
calc_aedes_recruitment_hill <- function(
  rain_14d,
  r_threshold = 50,
  r_peak = 130,
  r_flush = 200,
  base_frac = 0.05
) {
  rain_response <- ifelse(
    rain_14d < r_threshold,
    0,
    ifelse(
      rain_14d <= r_peak,
      (rain_14d - r_threshold) / (r_peak - r_threshold),
      pmax(0, (r_flush - rain_14d) / (r_flush - r_peak))
    )
  )
  base_frac + (1 - base_frac) * rain_response
}

#' Aedes recruitment scaling (0 to 1): smooth Briere-like response
calc_aedes_recruitment_briere <- function(
  rain_14d,
  r_threshold = 50,
  r_flush = 200,
  base_frac = 0.05
) {
  support_width <- max(r_flush - r_threshold, .Machine$double.eps)
  rain_response <- ifelse(
    rain_14d <= r_threshold | rain_14d >= r_flush,
    0,
    (rain_14d - r_threshold) * sqrt(pmax(r_flush - rain_14d, 0))
  )
  max_response <- 2 * support_width^(3 / 2) / (3 * sqrt(3))
  rain_response <- rain_response / max_response
  base_frac + (1 - base_frac) * pmax(0, pmin(1, rain_response))
}

#' Dispatch to the selected Aedes rainfall-response candidate
calc_aedes_recruitment <- function(
  rain_14d,
  candidate = AEDES_CURVE_CANDIDATE,
  r_flush = 200,
  base_frac = 0.05
) {
  candidate <- match.arg(candidate, AEDES_CURVE_CANDIDATES)

  if (candidate == "baseline_hill") {
    return(calc_aedes_recruitment_hill(
      rain_14d = rain_14d,
      r_threshold = 50,
      r_peak = 130,
      r_flush = r_flush,
      base_frac = base_frac
    ))
  }

  if (candidate == "early_flush_hill") {
    return(calc_aedes_recruitment_hill(
      rain_14d = rain_14d,
      r_threshold = 50,
      r_peak = 110,
      r_flush = r_flush,
      base_frac = base_frac
    ))
  }

  calc_aedes_recruitment_briere(
    rain_14d = rain_14d,
    r_threshold = 50,
    r_flush = r_flush,
    base_frac = base_frac
  )
}

#' Culex recruitment scaling (0 to 1)
#'
#' Hump-shaped rainfall response. Culex breed in permanent
#' and semi-permanent water (irrigation, livestock watering, lakes).
#' base_frac controls dry-season minimum. Decline above w_max due to flushing.
calc_culex_recruitment <- function(
  rain_14d,
  w_min = 5,
  w_max = 120,
  base_frac = 0.40
) {
  rain_response <- ifelse(
    rain_14d < w_min | rain_14d > w_max,
    0,
    4 * (rain_14d - w_min) * (w_max - rain_14d) / (w_max - w_min)^2
  )
  base_frac + (1 - base_frac) * pmax(0, rain_response)
}

# ==============================================================================
# Thermal parameters - Fischer et al. (2013)
# ==============================================================================

# Per-bite transmission differs by genus:
#   Aedes: highly competent, drives wet-season outbreaks
#   Culex: less efficient per bite (0.060), but persists year-round via abundance
thermal_params_aedes <- list(
  lf_a0 = 25.8,
  lf_a1 = 0.45, # Fischer (2013): lifespan = lf_a0 - lf_a1*T
  b_min = 9.6,
  b_slope = 0.0173, # Fischer (2013): biting rate = b_slope*(T - b_min)
  eip_max = 18.9,
  eip_slope = 0.30, # Fischer (2013): EIP = eip_max - eip_slope*T
  beta_hv = 0.65, # P(host infects Aedes | bite): very high (Turell et al. 2008)
  beta_vh = 0.50 # P(Aedes infects host | bite): very high
)

# Culex: abundant but less efficient per bite. Maintains dry-season persistence.
thermal_params_culex <- list(
  lf_a0 = 25.8,
  lf_a1 = 0.45,
  b_min = 9.6,
  b_slope = 0.0173,
  eip_max = 18.9,
  eip_slope = 0.30,
  beta_hv = 0.30, # P(host infects Culex | bite): moderate
  beta_vh = 0.20 # P(Culex infects host | bite): moderate
)

# ==============================================================================
# Temperature dependent parameter functions
# ==============================================================================

calc_lifespan <- function(Temp, species = "culex") {
  p <- if (species == "culex") thermal_params_culex else thermal_params_aedes
  pmax(1, p$lf_a0 - p$lf_a1 * Temp)
}
# Mortality rate is the inverse of lifespan.
calc_mortality_rate <- function(Temp, species = "culex") {
  1 / calc_lifespan(Temp, species)
}

calc_biting_rate <- function(Temp, species = "culex") {
  p <- if (species == "culex") thermal_params_culex else thermal_params_aedes
  pmax(0, p$b_slope * (Temp - p$b_min))
}
# EIP is the time until an infected mosquito becomes infectious.
calc_eip <- function(Temp, species = "culex") {
  p <- if (species == "culex") thermal_params_culex else thermal_params_aedes
  pmax(1, p$eip_max - p$eip_slope * Temp)
}
# EIP rate is the inverse of EIP duration.
calc_eip_rate <- function(Temp, species = "culex") {
  1 / calc_eip(Temp, species)
}

get_transmission_prob <- function(species = "culex", direction = "hv") {
  p <- if (species == "culex") thermal_params_culex else thermal_params_aedes
  if (direction == "hv") p$beta_hv else p$beta_vh
}

# ==============================================================================
# Static (non-climate) parameters
# ==============================================================================

params_static <- list(
  # Host demography
  mu_h = 1 / (5 * 365), # Host mortality: ~5 year lifespan (1/day)

  # Host infection dynamics (from Fischer Table S1)
  sigma_h = 1 / 4, # Latent rate: ~4 day latent period (1/day)
  gamma_h = 1 / 6, # Recovery rate: ~6 day infectious period (1/day)
  alpha_h = 0.0, # Disease-induced mortality (not modelled)
  omega_h = 1 / (6 * 365), # Waning immunity: R->S, ~6 year half-life (Essa's 2026's paper))

  # Vector parameters (static fallback, used when temperature dependence is off)
  mu_A = 1 / 15, # Aedes mortality: 15-day lifespan
  sigma_A = 1 / 10, # Aedes EIP rate: 10-day EIP
  a_A = 1 / 4, # Aedes biting: every 4 days

  mu_C = 1 / 15, # Culex mortality: 15-day lifespan
  sigma_C = 1 / 10, # Culex EIP rate: 10-day EIP
  a_C = 1 / 4, # Culex biting: every 4 days

  # Transmission probabilities (from thermal params)
  beta_vh_A = thermal_params_aedes$beta_vh, # Aedes -> host
  beta_vh_C = thermal_params_culex$beta_vh, # Culex -> host
  beta_hv_A = thermal_params_aedes$beta_hv, # Host -> Aedes
  beta_hv_C = thermal_params_culex$beta_hv, # Host -> Culex

  # Local coupling strength: fraction of time cattle spend near ward boundaries
  # (shared grazing, watering points).
  eps_local = 0.15,

  # Annualised exogenous infection pressure eta_i.
  # This is a phenomenological background term that keeps low-level infection
  # pressure entering from outside the modelled system when the observed ward
  # movement network is incomplete. It is converted to a daily quantity
  # internally and added to effective host infectious prevalence seen by
  # vectors, rather than representing explicit cattle import events.
  intro_rate = 0.001, # default: 0.1% annual probability per host

  # Which thermal parameter set to use
  temp_species_A = "aedes",
  temp_species_C = "culex"
)

# Recruitment scaling: call the species-specific function directly.
# Each function returns 0-1 with base_frac built in.
calc_recruit_scale <- function(
  rain_14d,
  species = c("aedes", "culex"),
  aedes_candidate = AEDES_CURVE_CANDIDATE
) {
  species <- match.arg(species)
  if (species == "aedes") {
    calc_aedes_recruitment(rain_14d, candidate = aedes_candidate)
  } else {
    calc_culex_recruitment(rain_14d)
  }
}

# ------------------------------------------------------------------------------
# Maximum vectors per host (carrying capacity)
# ------------------------------------------------------------------------------
#' District-level carrying capacity.
#' Kilimanjaro (Hai, Moshi): higher rainfall (supports more vectors per host
#' throughout the year)
#' Babati (Manyara): drier, seasonal flooding (according to Jess), lower
#' baseline vector abundance but higher Aedes pulses during floods.
#'
#' Species differences come from recruitment dynamics (Aedes boom-bust vs
#' Culex persistent) and transmission competence (Aedes >> Culex per bite).
#' From 4-parameter grid search: district-level VPH, separate for Aedes and Culex.
VPH_A_KILI <- 6 # Aedes in Hai + Moshi
VPH_C_KILI <- 10 # Culex in Hai + Moshi
VPH_A_BABATI <- 8 # Aedes in Babati
VPH_C_BABATI <- 8 # Culex in Babati


# ==============================================================================
# Plotting functions
# ==============================================================================

#' Plot rainfall-dependent recruitment functions
plot_rainfall_curves <- function(climate_df = NULL, sample_ward = NULL) {
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  require(patchwork)

  rain_seq <- seq(0, 200, by = 1)
  func_df <- tibble(
    rain_14d = rain_seq,
    Aedes = calc_aedes_recruitment(rain_seq),
    Culex = calc_culex_recruitment(rain_seq)
  )

  p1_functions <- func_df |>
    pivot_longer(-rain_14d, names_to = "Species", values_to = "Recruitment") |>
    ggplot(aes(x = rain_14d, y = Recruitment, color = Species)) +
    geom_line(linewidth = 1.2) +
    geom_vline(
      xintercept = c(50, 130),
      linetype = "dashed",
      color = "firebrick",
      alpha = 0.5
    ) +
    geom_vline(
      xintercept = c(5, 60, 120),
      linetype = "dashed",
      color = "forestgreen",
      alpha = 0.5
    ) +
    scale_color_manual(values = cols_vec) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
    theme_bw() +
    labs(
      title = "A. Recruitment functions",
      x = "14-day accumulated rainfall (mm)",
      y = "Recruitment scaling (0-1)",
      color = NULL
    )

  if (is.null(climate_df)) {
    return(list(functions = p1_functions, combined = p1_functions))
  }

  if (is.null(sample_ward)) {
    sample_ward <- unique(climate_df$fullname)[1]
  }

  rain_ts <- climate_df |>
    filter(fullname == sample_ward) |>
    select(date, rain, rain_14d) |>
    filter(!is.na(date))

  p2_rainfall <- rain_ts |>
    pivot_longer(-date, names_to = "type", values_to = "mm") |>
    mutate(
      type = factor(
        type,
        levels = c("rain", "rain_14d"),
        labels = c("Daily", "14-day sum")
      )
    ) |>
    ggplot(aes(x = date, y = mm, color = type)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    scale_color_manual(
      values = c("Daily" = "steelblue", "14-day sum" = "darkblue")
    ) +
    theme_bw() +
    labs(
      title = "B. Example rainfall data",
      subtitle = paste("Ward:", sample_ward),
      x = NULL,
      y = "Rainfall (mm)",
      color = NULL
    )

  scaling_ts <- climate_df |>
    filter(fullname == sample_ward) |>
    select(date, scale_aedes, scale_culex) |>
    filter(!is.na(date))

  p3_scaling <- scaling_ts |>
    pivot_longer(-date, names_to = "species", values_to = "scaling") |>
    mutate(
      species = factor(
        species,
        levels = c("scale_aedes", "scale_culex"),
        labels = c("Aedes", "Culex")
      )
    ) |>
    ggplot(aes(x = date, y = scaling, color = species)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = cols_vec) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    labs(
      title = "C. Applied recruitment scaling",
      x = NULL,
      y = "Recruitment scaling (0-1)",
      color = NULL
    )

  combined <- p1_functions /
    p2_rainfall /
    p3_scaling +
    plot_annotation(title = "Rainfall-dependent mosquito recruitment")

  list(
    functions = p1_functions,
    rainfall = p2_rainfall,
    scaling = p3_scaling,
    combined = combined
  )
}

#' Plot all named Aedes rainfall-response candidates on the same axes
plot_aedes_curve_candidates <- function(
  rain_seq = seq(0, 300, by = 1),
  r_flush = 200,
  base_frac = 0.05
) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  candidate_df <- tidyr::expand_grid(
    rain_14d = rain_seq,
    candidate = AEDES_CURVE_CANDIDATES
  ) |>
    mutate(
      label = vapply(candidate, get_aedes_curve_label, character(1)),
      recruitment = mapply(
        FUN = function(rain_val, candidate_name) {
          calc_aedes_recruitment(
            rain_14d = rain_val,
            candidate = candidate_name,
            r_flush = r_flush,
            base_frac = base_frac
          )
        },
        rain_14d, candidate
      )
    )

  candidate_df |>
    ggplot(aes(x = rain_14d, y = recruitment, colour = label)) +
    geom_line(linewidth = 1.1) +
    theme_bw(base_size = 12) +
    labs(
      x = "14-day accumulated rainfall (mm)",
      y = "Recruitment scaling (0-1)",
      colour = NULL,
      title = "Aedes rainfall-response candidates"
    )
}

#' Plot thermal performance curves
plot_thermal_curves <- function(temp_range = seq(0, 35, by = 0.5)) {
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  require(patchwork)

  df <- tibble(Temp = temp_range) |>
    mutate(
      lf_culex = calc_lifespan(Temp, "culex"),
      lf_aedes = calc_lifespan(Temp, "aedes"),
      mu_culex = calc_mortality_rate(Temp, "culex"),
      mu_aedes = calc_mortality_rate(Temp, "aedes"),
      a_culex = calc_biting_rate(Temp, "culex"),
      a_aedes = calc_biting_rate(Temp, "aedes"),
      eip_culex = calc_eip(Temp, "culex"),
      eip_aedes = calc_eip(Temp, "aedes")
    )

  p1 <- df |>
    select(Temp, Culex = lf_culex, Aedes = lf_aedes) |>
    pivot_longer(-Temp, names_to = "species", values_to = "days") |>
    ggplot(aes(Temp, days, color = species)) +
    geom_line(linewidth = 1.2, show.legend = F) +
    scale_color_manual(values = cols_vec) +
    labs(title = "Lifespan", x = "Temperature (C)", y = "Days", color = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")

  p2 <- df |>
    select(Temp, Culex = mu_culex, Aedes = mu_aedes) |>
    pivot_longer(-Temp, names_to = "species", values_to = "rate") |>
    ggplot(aes(Temp, rate, color = species)) +
    geom_line(linewidth = 1.2, show.legend = F) +
    scale_color_manual(values = cols_vec) +
    labs(
      title = "Mortality rate",
      x = "Temperature (C)",
      y = "1/day",
      color = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  p3 <- df |>
    select(Temp, Rate = a_culex) |>
    ggplot(aes(Temp, Rate)) +
    geom_line(linewidth = 1.2, color = "darkgreen") +
    labs(
      title = "Biting rate",
      x = "Temperature (C)",
      y = "Bites/day"
    ) +
    theme_bw()

  p4 <- df |>
    select(Temp, Culex = eip_culex, Aedes = eip_aedes) |>
    pivot_longer(-Temp, names_to = "species", values_to = "days") |>
    ggplot(aes(Temp, days, color = species)) +
    geom_line(linewidth = 1.2, show.legend = F) +
    scale_color_manual(values = cols_vec) +
    labs(
      title = "Extrinsic Incubation Period (EIP)",
      x = "Temperature (C)",
      y = "Days",
      color = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  (p1 + p2) /
    (p3 + p4) +
    plot_annotation(title = "Temperature-dependent vector parameters")
}


# ==============================================================================
# Climate overview plotting function
# ==============================================================================

#' Plot annual climate overview for hindcast validation
#'
#' Creates a 3-panel figure showing yearly rainfall, temperature, and NDVI
#' anomalies by district. Vertical reference lines mark known outbreak years.
plot_climate_overview <- function(
  ward_climate,
  ward_district,
  ndvi_df = NULL,
  outbreak_years = c(2006, 2007)
) {
  require(ggplot2)
  require(dplyr)
  require(patchwork)

  # Annual rainfall totals by district
  annual_rain <- ward_climate |>
    mutate(year = year(date), district = ward_district[fullname]) |>
    filter(!is.na(district)) |>
    group_by(district, fullname, year) |>
    summarise(annual_rain_mm = sum(rain, na.rm = TRUE), .groups = "drop")

  rain_district <- annual_rain |>
    group_by(district, year) |>
    summarise(
      mean_rain = mean(annual_rain_mm),
      lo = quantile(annual_rain_mm, 0.25),
      hi = quantile(annual_rain_mm, 0.75),
      .groups = "drop"
    )

  p_rain <- rain_district |>
    ggplot(aes(x = year, y = mean_rain, color = district, fill = district)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    geom_vline(
      xintercept = outbreak_years,
      linetype = "dashed",
      color = "red",
      alpha = 0.6
    ) +
    scale_color_manual(values = cols_district) +
    scale_fill_manual(values = cols_district) +
    theme_bw() +
    labs(
      title = "A. Annual rainfall",
      x = NULL,
      y = "Total rainfall (mm/year)",
      color = "District",
      fill = "District"
    )

  # Annual mean temperature by district
  annual_temp <- ward_climate |>
    mutate(year = year(date), district = ward_district[fullname]) |>
    filter(!is.na(district)) |>
    group_by(district, fullname, year) |>
    summarise(mean_temp = mean(temp, na.rm = TRUE), .groups = "drop")

  temp_district <- annual_temp |>
    group_by(district, year) |>
    summarise(
      mean_temp = mean(mean_temp),
      lo = quantile(mean_temp, 0.25),
      hi = quantile(mean_temp, 0.75),
      .groups = "drop"
    )

  p_temp <- temp_district |>
    ggplot(aes(x = year, y = mean_temp, color = district, fill = district)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    geom_vline(
      xintercept = outbreak_years,
      linetype = "dashed",
      color = "red",
      alpha = 0.6
    ) +
    scale_color_manual(values = cols_district) +
    scale_fill_manual(values = cols_district) +
    theme_bw() +
    labs(
      title = "B. Annual mean temperature",
      x = NULL,
      y = "Temperature (\u00b0C)",
      color = "District",
      fill = "District"
    )

  # NDVI panel (if data available)
  if (!is.null(ndvi_df) && nrow(ndvi_df) > 0) {
    annual_ndvi <- ndvi_df |>
      mutate(year = year(date), district = ward_district[fullname]) |>
      filter(!is.na(district)) |>
      group_by(district, fullname, year) |>
      summarise(mean_ndvi = mean(ndvi, na.rm = TRUE), .groups = "drop")

    ndvi_district <- annual_ndvi |>
      group_by(district, year) |>
      summarise(
        mean_ndvi = mean(mean_ndvi),
        lo = quantile(mean_ndvi, 0.25),
        hi = quantile(mean_ndvi, 0.75),
        .groups = "drop"
      )

    p_ndvi <- ndvi_district |>
      ggplot(aes(x = year, y = mean_ndvi, color = district, fill = district)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.5) +
      geom_vline(
        xintercept = outbreak_years,
        linetype = "dashed",
        color = "red",
        alpha = 0.6
      ) +
      scale_color_manual(values = cols_district) +
      scale_fill_manual(values = cols_district) +
      theme_bw() +
      labs(
        title = "C. Annual mean NDVI",
        x = "Year",
        y = "NDVI",
        color = "District",
        fill = "District"
      )

    combined <- p_rain /
      p_temp /
      p_ndvi +
      plot_annotation(
        title = "Climate overview: 2000\u20132024"
      )
  } else {
    combined <- p_rain /
      p_temp +
      plot_annotation(title = "Climate overview: 2000\u20132024")
  }

  combined
}
