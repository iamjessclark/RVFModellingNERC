# ==============================================================================
# 03_helpers.R
#
# Post-processing and metric extraction functions for model output, plus the
# shared calibration-target and objective machinery used by both the cluster
# calibration scripts and the ABC rejection step. Separated from 02_model.R
# for clarity.
#
# Output extraction:
#   extract_comp_mat()              — extract compartment as [time x wards] matrix
#   calc_global_metrics_from_out()  — peak prevalence and cumulative incidence
#   build_district_ts()             — aggregate ward output to district level
#   calc_ward_endemic_prevalence()  — ward-level prevalence from state vector
#   extract_regional_ts()           — regional time series summary
#   extract_ward_seroprevalence()   — ward seroprevalence at specific time
#   summarise_annual_foi_from_Hh()  — annual FOI from per-genus cumulative hazard
#   get_ward_coords()               — ward coordinates for spatial plots
#
# Calibration scaffolding:
#   build_calibration_targets()     — observed targets shared across scripts
#   score_calibration_run()         — one-sample objective + diagnostics
#   relative_mse()                  — robust mean relative squared error
#   weighted_mean_safe()            — weighted mean tolerating NAs/empty input
#   safe_metric_name()              — CSV-friendly column suffixes
#   named_metric_values()           — helper for ward-level diagnostic columns
#   clean_ward_id()                 — lower-case + trim ward identifiers
# ==============================================================================

# ==============================================================================
# Metric extraction functions
# ==============================================================================

#' Extract compartment values as matrix [n_time x n_wards]
#'
#' @param out_df Output from run_ode
#' @param prefix Compartment prefix (e.g., "Sh", "IvA")
#' @param ward_ids Ward identifiers
#' @return Matrix [n_time x n_wards]
extract_comp_mat <- function(out_df, prefix, ward_ids) {
  cols <- nm(prefix, ward_ids)
  as.matrix(out_df[, cols, drop = FALSE])
}


#' Calculate global outbreak metrics
#'
#' @param out_df Output from run_ode
#' @param ward_ids Ward identifiers
#' @return List with peak_prev (maximum regional prevalence) and
#'         cum_incidence (total infections over simulation)
calc_global_metrics_from_out <- function(out_df, ward_ids) {
  Sh <- extract_comp_mat(out_df, "Sh", ward_ids)
  Eh <- extract_comp_mat(out_df, "Eh", ward_ids)
  Ih <- extract_comp_mat(out_df, "Ih", ward_ids)
  Rh <- extract_comp_mat(out_df, "Rh", ward_ids)
  Ch <- extract_comp_mat(out_df, "Ch", ward_ids)

  Nh_tot <- rowSums(Sh + Eh + Ih + Rh)
  Ih_tot <- rowSums(Ih)

  prevalence_ts <- Ih_tot / pmax(Nh_tot, 1)
  peak_prev <- max(prevalence_ts, na.rm = TRUE)
  cum_incidence <- sum(Ch[nrow(Ch), ], na.rm = TRUE)

  list(peak_prev = peak_prev, cum_incidence = cum_incidence)
}


#' Build district-level time series
#'
#' Aggregates ward-level outputs to district level
#'
#' @param out_df Output from run_ode
#' @param dates Date sequence
#' @param ward_ids Ward identifiers
#' @param ward_district Named vector mapping ward to district
#' @return Long-format data frame with district-level metrics
build_district_ts <- function(out_df, dates, ward_ids, ward_district) {
  # Create district aggregation matrix
  district_fac <- factor(ward_district[ward_ids])
  D <- model.matrix(~ 0 + district_fac)
  colnames(D) <- gsub("^district_fac", "", colnames(D))

  # Aggregate compartments to district level
  Sh <- extract_comp_mat(out_df, "Sh", ward_ids) %*% D
  Eh <- extract_comp_mat(out_df, "Eh", ward_ids) %*% D
  Ih <- extract_comp_mat(out_df, "Ih", ward_ids) %*% D
  Rh <- extract_comp_mat(out_df, "Rh", ward_ids) %*% D
  Ch <- extract_comp_mat(out_df, "Ch", ward_ids) %*% D

  # Calculate derived metrics
  Nh <- Sh + Eh + Ih + Rh
  prevalence <- Ih / pmax(Nh, 1)
  seroprev <- Rh / pmax(Nh, 1)
  daily_inc <- apply(Ch, 2, function(x) pmax(0, c(x[1], diff(x))))

  # Combine into long-format data frame
  tibble(date = dates, time = seq_along(dates)) |>
    bind_cols(
      as_tibble(prevalence) |> setNames(paste0("prev__", colnames(prevalence))),
      as_tibble(seroprev) |> setNames(paste0("sero__", colnames(seroprev))),
      as_tibble(Ch) |> setNames(paste0("cum__", colnames(Ch))),
      as_tibble(daily_inc) |> setNames(paste0("inc__", colnames(daily_inc)))
    ) |>
    pivot_longer(
      -c(date, time),
      names_to = c("metric", "district"),
      names_sep = "__",
      values_to = "value"
    ) |>
    pivot_wider(names_from = metric, values_from = value)
}


#' Calculate endemic prevalence by ward from final state
#'
#' @param state State vector (e.g., from end of burn-in)
#' @param ward_ids Ward identifiers
#' @return Data frame with ward-level prevalence and seroprevalence
calc_ward_endemic_prevalence <- function(state, ward_ids) {
  tibble(ward = ward_ids) |>
    rowwise() |>
    mutate(
      Sh = state[paste0("Sh__", ward)],
      Eh = state[paste0("Eh__", ward)],
      Ih = state[paste0("Ih__", ward)],
      Rh = state[paste0("Rh__", ward)],
      Nh = max(Sh + Eh + Ih + Rh, 1),
      prevalence = Ih / Nh,
      seroprevalence = Rh / Nh
    ) |>
    ungroup()
}


# ==============================================================================
# Annual FOI from cumulative hazard
# ==============================================================================

#' Calculate annual FOI probability from cumulative hazard
#'
#' The ODE tracks cumulative host hazard split by vector genus
#' (Hh_A for Aedes-mediated, Hh_C for Culex-mediated). Total host hazard
#' is Hh_A + Hh_C, and annual infection probability is 1 - exp(-hazard).
#'
#' Returns the total hazard and probability (backward-compatible columns
#' annual_foi_hazard and annual_foi_prob), plus per-genus contributions
#' (annual_foi_hazard_A/_C and annual_foi_prob_A/_C). The per-genus
#' probabilities are computed AS IF that genus acted alone, i.e.
#' 1 - exp(-Hh_g). They are not the genus' share of the total
#' probability; for that, use the hazard scale where shares are additive.
#'
#' @param out_df Output from run_ode (must have date column)
#' @param ward_ids Ward identifiers
#' @param date_start Start of period
#' @param date_end End of period
#' @return Data frame with ward, total + per-genus hazard and infection probability
summarise_annual_foi_from_Hh <- function(out_df, ward_ids, date_start, date_end) {
  stopifnot("date" %in% names(out_df))

  i_start <- which(out_df$date == date_start)
  i_end   <- which(out_df$date == date_end)

  if (length(i_start) != 1 || length(i_end) != 1) {
    stop("Dates not found in output: ", date_start, " to ", date_end)
  }

  # Read genus-specific cumulative hazard. Total = Hh_A + Hh_C.
  Hh_A_mat <- extract_comp_mat(out_df, "Hh_A", ward_ids)
  Hh_C_mat <- extract_comp_mat(out_df, "Hh_C", ward_ids)

  H_A <- as.numeric(Hh_A_mat[i_end, ] - Hh_A_mat[i_start, ])
  H_C <- as.numeric(Hh_C_mat[i_end, ] - Hh_C_mat[i_start, ])
  H_T <- H_A + H_C

  # n_days kept for downstream code that reports the period length; the
  # off-by-one is a minor lsoda-output convention noted in the original
  # implementation.
  n_days <- as.integer(date_end - date_start) + 1L

  tibble(
    ward = ward_ids,
    date_start = date_start,
    date_end   = date_end,
    n_days     = n_days,
    # Backward-compatible total columns (used everywhere downstream).
    annual_foi_hazard = H_T,
    annual_foi_prob   = 1 - exp(-H_T),
    # New per-genus columns. Hazard is additive across genera; probability
    # here is computed as if each genus acted alone (genus-only counterfactual
    # equivalent on the probability scale).
    annual_foi_hazard_A = H_A,
    annual_foi_hazard_C = H_C,
    annual_foi_prob_A   = 1 - exp(-H_A),
    annual_foi_prob_C   = 1 - exp(-H_C)
  )
}


#' Summarise interval seroconversion in a passive baseline-seronegative cohort
#'
#' The follow-up survey counts new seropositive animals among the longitudinal
#' cohort, so the closest model-equivalent estimand is not cumulative Rh/Nh at
#' follow-up. This helper starts an observation cohort at `date_start` among
#' animals that are not in R, then lets that cohort progress through the host
#' infection and recovery stages using the model's daily host FOI. Cohort
#' members entering the absorbing "new R" class are counted as seroconverted.
#'
#' This is an observation model only: it does not feed back into transmission,
#' does not add births, and does not count animals that were already R at
#' baseline.
#'
#' @param out_df Output from run_ode (must have date column and Hh_A/Hh_C)
#' @param ward_ids Ward identifiers
#' @param date_start Start of follow-up window
#' @param date_end End of follow-up window
#' @param params Parameter list containing sigma_h and gamma_h
#' @param n_substeps Number of within-day substeps for the passive cohort
#' @param initial_state How to initialise the baseline-seronegative cohort.
#'   "nonR_distribution" uses the model S/E/I mix among non-R animals at
#'   baseline; "susceptible_only" starts all baseline-seronegative animals in S.
#' @return Data frame with ward-level seroconversion probability
summarise_interval_seroconversion_from_R_flow <- function(
    out_df, ward_ids, date_start, date_end, params,
    n_substeps = 4,
    initial_state = c("nonR_distribution", "susceptible_only")) {
  stopifnot("date" %in% names(out_df))
  initial_state <- match.arg(initial_state)

  i_start <- which(out_df$date == date_start)
  i_end   <- which(out_df$date == date_end)
  if (length(i_start) != 1 || length(i_end) != 1 || i_end <= i_start) {
    stop("Dates not found, or invalid interval: ", date_start, " to ", date_end)
  }
  if (is.null(params$sigma_h) || is.null(params$gamma_h)) {
    stop("params must contain sigma_h and gamma_h.")
  }

  Sh <- extract_comp_mat(out_df, "Sh", ward_ids)
  Eh <- extract_comp_mat(out_df, "Eh", ward_ids)
  Ih <- extract_comp_mat(out_df, "Ih", ward_ids)
  Rh <- extract_comp_mat(out_df, "Rh", ward_ids)
  Hh_A <- extract_comp_mat(out_df, "Hh_A", ward_ids)
  Hh_C <- extract_comp_mat(out_df, "Hh_C", ward_ids)

  H_total <- Hh_A + Hh_C
  H_window <- as.numeric(H_total[i_end, ] - H_total[i_start, ])
  daily_hazard <- pmax(
    H_total[(i_start + 1):i_end, , drop = FALSE] -
      H_total[i_start:(i_end - 1), , drop = FALSE],
    0
  )

  nonR_start <- as.numeric(Sh[i_start, ] + Eh[i_start, ] + Ih[i_start, ])
  if (initial_state == "susceptible_only") {
    S <- rep(1, length(ward_ids))
    E <- I <- rep(0, length(ward_ids))
  } else {
    denom <- pmax(nonR_start, .Machine$double.eps)
    S <- as.numeric(Sh[i_start, ]) / denom
    E <- as.numeric(Eh[i_start, ]) / denom
    I <- as.numeric(Ih[i_start, ]) / denom
  }
  R_new <- rep(0, length(ward_ids))

  n_substeps <- max(1L, as.integer(n_substeps))
  dt <- 1 / n_substeps
  p_EI <- 1 - exp(-params$sigma_h * dt)
  p_IR <- 1 - exp(-params$gamma_h * dt)

  for (day_idx in seq_len(nrow(daily_hazard))) {
    lambda_day <- as.numeric(daily_hazard[day_idx, ])
    for (sub_idx in seq_len(n_substeps)) {
      p_SE <- 1 - exp(-(lambda_day * dt))
      new_E <- p_SE * S
      new_I <- p_EI * E
      new_R <- p_IR * I

      S <- pmax(S - new_E, 0)
      E <- pmax(E + new_E - new_I, 0)
      I <- pmax(I + new_I - new_R, 0)
      R_new <- pmin(R_new + pmax(new_R, 0), 1)
    }
  }

  tibble(
    ward = ward_ids,
    date_start = date_start,
    date_end = date_end,
    n_days = as.integer(date_end - date_start) + 1L,
    nonR_start = nonR_start,
    annual_foi_hazard = H_window,
    annual_foi_prob = 1 - exp(-H_window),
    seroconv_prob_newR = pmin(pmax(R_new, 0), 1),
    observation_model = paste0("new_R_", initial_state)
  )
}


# ==============================================================================
# Coordinate extraction (for spatial plots)
# ==============================================================================

#' Pick first matching column from candidates (case-insensitive)
pick_first_col <- function(df, candidates) {
  nms <- names(df)
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    return(NA_character_)
  }
  nms[idx[1]]
}


#' Get ward coordinates from ward_info data frame
#'
#' @param ward_info Data frame with ward info including coordinates
#' @return Data frame with ward, lon, lat columns; or NULL if not available
get_ward_coords <- function(ward_info) {
  lon_col <- pick_first_col(ward_info, c("Longitude", "longitude", "lon", "x"))
  lat_col <- pick_first_col(ward_info, c("Latitude", "latitude", "lat", "y"))
  id_col <- pick_first_col(ward_info, c("clean_name", "fullname"))

  if (is.na(lon_col) || is.na(lat_col)) {
    return(NULL)
  }

  ward_info |>
    mutate(
      ward = if (!is.na(id_col)) .data[[id_col]] else .data[["fullname"]],
      ward = clean_id(ward)
    ) |>
    transmute(ward, lon = .data[[lon_col]], lat = .data[[lat_col]]) |>
    distinct()
}


# ==============================================================================
# Regional time series extraction
# ==============================================================================

#' Extract region-wide time series from ODE output
#'
#' Aggregates all wards into a single regional summary with host and vector
#' compartments, prevalence, seroprevalence, and vectors-per-host ratios.
#'
#' @param out_df Output from run_ode (data frame with state columns and date)
#' @param scenario_name Label for this scenario (e.g., "hindcast", "burn-in")
#' @return Data frame with regional totals and derived metrics
extract_regional_ts <- function(out_df, scenario_name) {
  out_df |>
    select(
      time, date,
      starts_with("Sh__"), starts_with("Eh__"),
      starts_with("Ih__"), starts_with("Rh__"),
      starts_with("SvA__"), starts_with("EvA__"), starts_with("IvA__"),
      starts_with("SvC__"), starts_with("EvC__"), starts_with("IvC__")
    ) |>
    pivot_longer(-c(time, date), names_to = c("compartment", "ward"), names_sep = "__") |>
    group_by(time, date, compartment) |>
    summarise(total = sum(value), .groups = "drop") |>
    pivot_wider(names_from = compartment, values_from = total) |>
    mutate(
      Nh = Sh + Eh + Ih + Rh,
      prevalence = Ih / pmax(Nh, 1),
      seroprevalence = Rh / pmax(Nh, 1),
      NvA = SvA + EvA + IvA,
      NvC = SvC + EvC + IvC,
      Aedes = NvA / pmax(Nh, 1),
      Culex = NvC / pmax(Nh, 1),
      scenario = scenario_name
    )
}


# ==============================================================================
# Survey data extraction (for calibration)
# ==============================================================================

#' Standardise ward identifiers outside the main data-prep pipeline
clean_ward_id <- function(x) tolower(trimws(x))


#' Make stable, CSV-friendly suffixes for saved calibration summary columns
safe_metric_name <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  y <- gsub("^_+|_+$", "", y)
  tolower(y)
}


#' Weighted mean that returns NA rather than failing on empty/missing inputs
weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w = w[ok])
}


#' Mean relative squared error with finite-value guarding
relative_mse <- function(model, observed) {
  err <- ((model - observed) / observed)^2
  err <- err[is.finite(err)]
  if (length(err) == 0) {
    return(Inf)
  }
  mean(err)
}


#' Return a list of named metric columns for survey-ward diagnostics
named_metric_values <- function(df, id_col, value_col, prefix, ids = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    if (is.null(ids)) {
      return(list())
    }
    out <- as.list(rep(NA_real_, length(ids)))
    names(out) <- paste0(prefix, safe_metric_name(ids))
    return(out)
  }

  values <- setNames(df[[value_col]], df[[id_col]])
  if (is.null(ids)) {
    ids <- names(values)
  }
  out <- as.list(as.numeric(values[ids]))
  names(out) <- paste0(prefix, safe_metric_name(ids))
  out
}


#' Clamp probabilities away from exact 0/1 for likelihood calculations
clamp_probability <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}


#' Binomial deviance against the saturated model
binom_deviance_safe <- function(k, n, p) {
  if (any(!is.finite(p))) {
    stop("Non-finite model probability supplied to binomial deviance.")
  }
  p <- clamp_probability(p)
  obs_p <- ifelse(n > 0, k / n, NA_real_)
  nll <- -sum(dbinom(k, size = n, prob = p, log = TRUE), na.rm = TRUE)
  sat <- -sum(dbinom(k, size = n, prob = obs_p, log = TRUE), na.rm = TRUE)
  2 * (nll - sat)
}


#' Binomial negative log-likelihood
binom_nll_safe <- function(k, n, p) {
  if (any(!is.finite(p))) {
    stop("Non-finite model probability supplied to binomial likelihood.")
  }
  -sum(dbinom(k, size = n, prob = clamp_probability(p), log = TRUE), na.rm = TRUE)
}


#' Beta-binomial negative log-likelihood with fixed precision
betabinom_nll_safe <- function(k, n, p, phi = 50) {
  if (any(!is.finite(p))) {
    stop("Non-finite model probability supplied to beta-binomial likelihood.")
  }
  p <- clamp_probability(p)
  alpha <- p * phi
  beta <- (1 - p) * phi
  -sum(lchoose(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(alpha, beta),
       na.rm = TRUE)
}


#' Build direct survey-count targets for likelihood-based calibration
#'
#' These are the raw survey estimands, not Paul/Jess smoothed FOI summaries:
#'   - 2023 baseline:  n_pos_2023 / n_tested_2023
#'   - 2024 follow-up: n_seroconv_2024 / n_resampled_2024
#'
#' @param survey_path Path to survey validation RDS
#' @return List with survey rows for baseline serology and follow-up conversion
build_direct_count_targets <- function(survey_path = "ward_info/survey_validation.rds") {
  survey_data <- readRDS(survey_path)
  ws <- survey_data$ward_sero |>
    mutate(
      clean_name = clean_ward_id(clean_name),
      district = sapply(strsplit(clean_name, "/"), function(x) x[2])
    )

  survey_2023 <- ws |>
    filter(!is.na(n_tested_2023), n_tested_2023 > 0, !is.na(n_pos_2023))

  survey_2024 <- ws |>
    filter(!is.na(n_resampled_2024), n_resampled_2024 > 0,
           !is.na(n_seroconv_2024))

  if (nrow(survey_2023) == 0 || nrow(survey_2024) == 0) {
    stop("Direct survey-count targets are empty; check survey_validation.rds.")
  }

  district_2023 <- survey_2023 |>
    group_by(district) |>
    summarise(
      obs = sum(n_pos_2023, na.rm = TRUE) /
        sum(n_tested_2023, na.rm = TRUE),
      n = sum(n_tested_2023, na.rm = TRUE),
      k = sum(n_pos_2023, na.rm = TRUE),
      n_wards = n(),
      .groups = "drop"
    ) |>
    arrange(desc(obs))

  district_2024 <- survey_2024 |>
    group_by(district) |>
    summarise(
      obs = sum(n_seroconv_2024, na.rm = TRUE) /
        sum(n_resampled_2024, na.rm = TRUE),
      n = sum(n_resampled_2024, na.rm = TRUE),
      k = sum(n_seroconv_2024, na.rm = TRUE),
      n_wards = n(),
      .groups = "drop"
    ) |>
    arrange(desc(obs))

  list(
    ws = ws,
    survey_2023 = survey_2023,
    survey_2024 = survey_2024,
    wards_2023 = survey_2023$clean_name,
    wards_2024 = survey_2024$clean_name,
    district_2023 = district_2023,
    district_2024 = district_2024,
    obs_regional_sero23 =
      sum(survey_2023$n_pos_2023, na.rm = TRUE) /
      sum(survey_2023$n_tested_2023, na.rm = TRUE),
    obs_regional_conv24 =
      sum(survey_2024$n_seroconv_2024, na.rm = TRUE) /
      sum(survey_2024$n_resampled_2024, na.rm = TRUE),
    observed_rank_2024 = paste(district_2024$district, collapse = ">")
  )
}


#' Score one run against direct survey counts
#'
#' @param sero Ward-level seroprevalence from extract_ward_seroprevalence()
#' @param prob_2024_8mo Ward-level 2024 follow-up model-equivalent probability.
#'   By default this reads annual_foi_prob, but callers can set conv24_col to
#'   score a seroconversion-specific observation model such as
#'   seroconv_prob_newR.
#' @param direct_targets Output from build_direct_count_targets()
#' @param conv24_col Column in prob_2024_8mo to compare with observed 2024
#'   seroconversion counts
#' @return List of direct-count objectives and diagnostics
score_direct_count_run <- function(sero, prob_2024_8mo,
                                   direct_targets = build_direct_count_targets(),
                                   conv24_col = "annual_foi_prob") {
  survey_2023 <- direct_targets$survey_2023
  survey_2024 <- direct_targets$survey_2024

  if (!conv24_col %in% names(prob_2024_8mo)) {
    stop("Column not found in 2024 model probabilities: ", conv24_col)
  }

  sero_model <- sero |>
    transmute(ward = clean_ward_id(ward_id), model_sero23 = seroprev)

  conv_model <- prob_2024_8mo |>
    transmute(ward = clean_ward_id(ward), model_conv24 = .data[[conv24_col]])

  df23 <- survey_2023 |>
    left_join(sero_model, by = c("clean_name" = "ward"))

  df24 <- survey_2024 |>
    left_join(conv_model, by = c("clean_name" = "ward"))

  missing_2023 <- df23$clean_name[!is.finite(df23$model_sero23)]
  missing_2024 <- df24$clean_name[!is.finite(df24$model_conv24)]
  if (length(missing_2023) > 0 || length(missing_2024) > 0) {
    stop(
      "Missing/non-finite model probabilities for direct-count scoring: ",
      paste(c(missing_2023, missing_2024), collapse = ", ")
    )
  }

  dev23 <- binom_deviance_safe(
    df23$n_pos_2023, df23$n_tested_2023, df23$model_sero23
  )
  dev24 <- binom_deviance_safe(
    df24$n_seroconv_2024, df24$n_resampled_2024, df24$model_conv24
  )
  nll23 <- binom_nll_safe(
    df23$n_pos_2023, df23$n_tested_2023, df23$model_sero23
  )
  nll24 <- binom_nll_safe(
    df24$n_seroconv_2024, df24$n_resampled_2024, df24$model_conv24
  )
  bb23 <- betabinom_nll_safe(
    df23$n_pos_2023, df23$n_tested_2023, df23$model_sero23
  )
  bb24 <- betabinom_nll_safe(
    df24$n_seroconv_2024, df24$n_resampled_2024, df24$model_conv24
  )

  n23 <- sum(df23$n_tested_2023, na.rm = TRUE)
  n24 <- sum(df24$n_resampled_2024, na.rm = TRUE)

  district_2023 <- df23 |>
    group_by(district) |>
    summarise(
      obs = sum(n_pos_2023, na.rm = TRUE) / sum(n_tested_2023, na.rm = TRUE),
      model = weighted_mean_safe(model_sero23, n_tested_2023),
      n = sum(n_tested_2023, na.rm = TRUE),
      k = sum(n_pos_2023, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(model))

  district_2024 <- df24 |>
    group_by(district) |>
    summarise(
      obs = sum(n_seroconv_2024, na.rm = TRUE) /
        sum(n_resampled_2024, na.rm = TRUE),
      model = weighted_mean_safe(model_conv24, n_resampled_2024),
      n = sum(n_resampled_2024, na.rm = TRUE),
      k = sum(n_seroconv_2024, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(model))

  district_value <- function(df, district_name, value_col) {
    value <- df[[value_col]][df$district == district_name]
    if (length(value) == 0) {
      return(NA_real_)
    }
    value[[1]]
  }

  district_named_values <- function(df, prefix) {
    out <- c()
    for (district_name in sort(unique(df$district))) {
      out[paste0("direct_obs_", prefix, "_", district_name)] <-
        district_value(df, district_name, "obs")
      out[paste0("direct_model_", prefix, "_", district_name)] <-
        district_value(df, district_name, "model")
    }
    as.list(out)
  }

  rank_2023 <- paste(district_2023$district, collapse = ">")
  rank_2024 <- paste(district_2024$district, collapse = ">")

  b23 <- district_value(district_2023, "babati", "model")
  b23_obs <- district_value(district_2023, "babati", "obs")
  b24 <- district_value(district_2024, "babati", "model")

  ward_diagnostics <- c(
    named_metric_values(
      df23, "clean_name", "model_sero23",
      "ward_sero23_", direct_targets$wards_2023
    ),
    named_metric_values(
      df24, "clean_name", "model_conv24",
      "ward_seroconv_2024_8mo_", direct_targets$wards_2024
    )
  )

  c(
    list(
      obj_direct_counts = (dev23 / n23) + (dev24 / n24),
      direct_obj_equal_year = (dev23 / n23) + (dev24 / n24),
      direct_binom_dev_2023 = dev23,
      direct_binom_dev_2024 = dev24,
      direct_binom_dev_total = dev23 + dev24,
      direct_binom_dev_2023_per_animal = dev23 / n23,
      direct_binom_dev_2024_per_animal = dev24 / n24,
      direct_binom_nll_2023 = nll23,
      direct_binom_nll_2024 = nll24,
      direct_binom_nll_total = nll23 + nll24,
      direct_betabinom_phi50_nll_2023 = bb23,
      direct_betabinom_phi50_nll_2024 = bb24,
      direct_betabinom_phi50_nll_total = bb23 + bb24,
      direct_obs_sero23 = direct_targets$obs_regional_sero23,
      direct_model_sero23 = weighted_mean_safe(df23$model_sero23, df23$n_tested_2023),
      direct_obs_conv24 = direct_targets$obs_regional_conv24,
      direct_model_conv24 = weighted_mean_safe(df24$model_conv24, df24$n_resampled_2024),
      direct_rank_2023 = rank_2023,
      direct_rank_2024 = rank_2024,
      direct_observed_rank_2024 = direct_targets$observed_rank_2024,
      direct_babati_gain = b24 - b23,
      direct_plausible_low_babati_high_2024 =
        isTRUE(b23 <= b23_obs + 0.03) && isTRUE(b24 >= 0.10),
      direct_babati_top_2024 = startsWith(rank_2024, "babati"),
      direct_babati_high_2024 = isTRUE(b24 >= 0.10)
    ),
    district_named_values(district_2023, "sero23"),
    district_named_values(district_2024, "conv24"),
    ward_diagnostics
  )
}


#' Build calibration targets shared across grid, cluster, and ABC scripts
#'
#' @param survey_path Path to survey validation RDS
#' @param ward_summary_path Path to ward summary RDS
#' @param N_star Optional named vector of ward cattle populations
#' @param exclude_pattern Regex for wards to exclude from the study region
#' @return List with observed targets and surveyed wards
build_calibration_targets <- function(survey_path = "ward_info/survey_validation.rds",
                                      ward_summary_path = "ward_info/ward_summary_df.rds",
                                      N_star = NULL,
                                      exclude_pattern = "babati.*urban") {
  # If the caller has not supplied N_star, rebuild it exactly the same way as
  # the modelling pipeline: lower-case ward IDs, exclude Babati urban wards,
  # and floor cattle population at 1000 to match 01_data_prep.R.
  if (is.null(N_star)) {
    ward_info <- readRDS(ward_summary_path) |>
      mutate(clean_name = clean_ward_id(fullname)) |>
      filter(!grepl(exclude_pattern, clean_name)) |>
      mutate(N_star = pmax(Cattle_pop, 1000))

    N_star <- setNames(ward_info$N_star, ward_info$clean_name)
  }

  # Load the survey table used for calibration and standardise ward IDs so the
  # observed data join to model outputs and N_star lookups without silent mismatches.
  survey_data <- readRDS(survey_path)
  ws <- survey_data$ward_sero |>
    mutate(
      clean_name = clean_ward_id(clean_name),
      district = sapply(strsplit(clean_name, "/"), function(x) x[2])
    )

  # Fail early if any survey ward is not present in the N_star lookup. This is
  # safer than letting weighted means quietly return NA later on.
  missing_weights <- setdiff(unique(ws$clean_name), names(N_star))
  if (length(missing_weights) > 0) {
    stop("Missing N_star values for survey wards: ", paste(missing_weights, collapse = ", "))
  }

  # Restrict FOI-based comparisons to wards with observed 2023 FOI rate.
  surveyed_wards <- ws |>
    filter(!is.na(foi_rate_2023), foi_rate_2023 > 0) |>
    pull(clean_name)

  sero_wards <- ws |>
    filter(!is.na(seroprev_2023), !is.na(n_tested_2023), n_tested_2023 > 0) |>
    pull(clean_name)

  # Calibration is on probability scale.
  #   2023 target = 1 - exp(-foi_rate_2023)            (annual infection prob)
  #   2024 target = 1 - exp(-foi_rate_2024 * 8/12)     (8-month conversion prob)
  # The survey-prep already materialised these as seroconv_hat_{2023,2024_8mo}.
  # Compare to model `1 - exp(-H over matched window)` (`infection_prob` /
  # `annual_foi_prob` from summarise_annual_foi_from_Hh).
  #
  # Simple ward-mean (not cattle-weighted) is intentional. Paul's catalytic
  # model is fit at the village level with hierarchical RE + spatial field;
  # the village FOIs were already aggregated to ward by tested-cattle weights
  # in prepare_survey_data.R. Cattle-weighting the district summary was
  # tested as a sensitivity — qualitative conclusions unchanged.
  obs_district_seroconv_2023 <- ws |>
    filter(!is.na(seroconv_hat_2023), seroconv_hat_2023 > 0) |>
    group_by(district) |>
    summarise(obs = mean(seroconv_hat_2023), .groups = "drop") |>
    arrange(desc(obs))

  obs_district_seroconv_2024_8mo <- ws |>
    filter(!is.na(seroconv_hat_2024_8mo), seroconv_hat_2024_8mo > 0) |>
    group_by(district) |>
    summarise(obs = mean(seroconv_hat_2024_8mo), .groups = "drop") |>
    arrange(desc(obs))

  obs_regional_seroconv_2023      <- mean(ws$seroconv_hat_2023, na.rm = TRUE)
  obs_regional_seroconv_2024_8mo  <- mean(ws$seroconv_hat_2024_8mo, na.rm = TRUE)

  # Seroprevalence is population-weighted, not a simple mean, because the model
  # serology summary is also population-weighted across wards.
  obs_sero <- weighted.mean(
    ws$seroprev_2023,
    w = N_star[ws$clean_name],
    na.rm = TRUE
  )

  # Survey-matched serology uses the ward sample sizes as weights. This is kept
  # alongside the population-weighted regional serology so we can switch between
  # population and survey estimands without rerunning the ODE simulations.
  obs_sero_sample <- weighted_mean_safe(ws$seroprev_2023, ws$n_tested_2023)

  obs_district_sero <- ws |>
    filter(clean_name %in% sero_wards) |>
    group_by(district) |>
    summarise(
      obs_sero = weighted_mean_safe(seroprev_2023, n_tested_2023),
      obs_sero_pop = weighted_mean_safe(seroprev_2023, N_star[clean_name]),
      n_tested = sum(n_tested_2023, na.rm = TRUE),
      n_wards = n(),
      .groups = "drop"
    ) |>
    arrange(district)

  # Save the observed district ordering explicitly. ABC can then compare rank
  # strings without rebuilding the ordering logic in a separate script.
  observed_rank      <- paste(obs_district_seroconv_2023$district,     collapse = ">")
  observed_rank_2024 <- paste(obs_district_seroconv_2024_8mo$district, collapse = ">")

  list(
    ws = ws,
    surveyed_wards = surveyed_wards,
    sero_wards     = sero_wards,
    obs_district_seroconv_2023     = obs_district_seroconv_2023,
    obs_district_seroconv_2024_8mo = obs_district_seroconv_2024_8mo,
    obs_district_sero              = obs_district_sero,
    obs_regional_seroconv_2023     = obs_regional_seroconv_2023,
    obs_regional_seroconv_2024_8mo = obs_regional_seroconv_2024_8mo,
    obs_sero        = obs_sero,
    obs_sero_pop    = obs_sero,
    obs_sero_sample = obs_sero_sample,
    observed_rank      = observed_rank,
    observed_rank_2024 = observed_rank_2024,
    target = c(
      seroconv_2023     = obs_regional_seroconv_2023,
      seroconv_2024_8mo = obs_regional_seroconv_2024_8mo,
      sero              = obs_sero,
      sero_sample       = obs_sero_sample,
      corr              = 1.0
    )
  )
}


#' Score one calibration run against the shared targets
#'
#' Both arguments are ward-level model summaries from
#' summarise_annual_foi_from_Hh() over the relevant calibration window. The
#' `annual_foi_prob` column on each is `1 - exp(-H)` over that window, i.e. an
#' **infection probability**, not a hazard rate. For v3:
#'   - prob_2023:     annual prob over ~2022-09-14 to 2023-09-14
#'   - prob_2024_8mo: 8-month prob over ~2023-09-27 to 2024-05-27
#' Targets are also probability-scale (see build_calibration_targets).
#'
#' @param prob_2023 Ward-level 12-month infection probability for 2023
#' @param sero Ward-level seroprevalence with district column
#' @param calibration_targets Output from build_calibration_targets()
#' @param N_star Named vector of cattle populations by ward
#' @param prob_2024_8mo Optional ward-level 8-month infection probability for 2024
#' @return List of model summaries, objective components, and total objective
score_calibration_run <- function(prob_2023, sero, calibration_targets, N_star,
                                  prob_2024_8mo = NULL) {
  surveyed_wards <- calibration_targets$surveyed_wards
  sero_wards     <- calibration_targets$sero_wards

  # Match model outputs to the ward subset used for the observed targets.
  prob_2023_surveyed     <- prob_2023 |> filter(ward %in% surveyed_wards)
  prob_2024_8mo_surveyed <- if (!is.null(prob_2024_8mo)) {
    prob_2024_8mo |> filter(ward %in% surveyed_wards)
  } else {
    NULL
  }
  sero_surveyed <- sero |> filter(ward_id %in% sero_wards)

  if (nrow(prob_2023_surveyed) == 0 || nrow(sero_surveyed) == 0) {
    stop("No surveyed wards found in model outputs for calibration scoring.")
  }

  # Regional model infection probability across surveyed wards.
  model_regional_seroconv_2023 <- mean(prob_2023_surveyed$annual_foi_prob)
  model_regional_seroconv_2024_8mo <- if (!is.null(prob_2024_8mo_surveyed) &&
                                         nrow(prob_2024_8mo_surveyed) > 0) {
    mean(prob_2024_8mo_surveyed$annual_foi_prob)
  } else {
    NA_real_
  }

  # Model seroprevalence uses the same cattle-population weights as the observed target.
  sero_weights <- N_star[sero_surveyed$ward_id]
  model_sero <- weighted_mean_safe(sero_surveyed$seroprev, sero_weights)

  sero_with_survey <- sero_surveyed |>
    left_join(
      calibration_targets$ws |>
        select(ward_id = clean_name, n_tested_2023),
      by = "ward_id"
    )
  model_sero_sample <- weighted_mean_safe(
    sero_with_survey$seroprev,
    sero_with_survey$n_tested_2023
  )

  # District means capture the broad spatial gradient we want the model to reproduce.
  # Simple ward-mean across surveyed wards, matching the observed target's
  # aggregation in build_calibration_targets() and Paul's reporting level.
  model_district_2023 <- prob_2023_surveyed |>
    group_by(district) |>
    summarise(prob = mean(annual_foi_prob), .groups = "drop") |>
    arrange(desc(prob))

  model_district_2024_8mo <- if (!is.null(prob_2024_8mo_surveyed) &&
                                 nrow(prob_2024_8mo_surveyed) > 0) {
    prob_2024_8mo_surveyed |>
      group_by(district) |>
      summarise(prob = mean(annual_foi_prob), .groups = "drop") |>
      arrange(desc(prob))
  } else {
    tibble(district = character(), prob = numeric())
  }

  model_district_sero <- sero_with_survey |>
    group_by(district) |>
    summarise(
      sero = weighted_mean_safe(seroprev, n_tested_2023),
      sero_pop = weighted_mean_safe(seroprev, N_star[ward_id]),
      .groups = "drop"
    ) |>
    arrange(district)

  # Ward-level correlation adds a finer-scale spatial-pattern check on top of
  # the district summaries and overall regional level. Compares model 2023
  # infection probability against the survey probability-scale target
  # (1 - exp(-foi_rate_2023)).
  matched <- inner_join(
    prob_2023_surveyed |> select(ward, model_prob = annual_foi_prob),
    calibration_targets$ws |>
      select(ward = clean_name, obs_prob = seroconv_hat_2023) |>
      filter(!is.na(obs_prob), obs_prob > 0),
    by = "ward"
  )
  corr <- if (nrow(matched) > 3) cor(matched$model_prob, matched$obs_prob) else 0

  district_err_df <- calibration_targets$obs_district_seroconv_2023 |>
    left_join(model_district_2023 |> rename(model_prob = prob), by = "district")

  district_sero_err_df <- calibration_targets$obs_district_sero |>
    select(district, obs_sero) |>
    left_join(model_district_sero |> rename(model_sero = sero), by = "district")

  # Objective components (relative-MSE on probability scale unless noted):
  #   foi_err       : overall 2023 infection probability mismatch
  #   dist_err      : 2023 district-level seroconv gradient mismatch
  #   foi_err_24    : overall 2024 8-month conversion mismatch (Inf if no 2024)
  #   dist_err_24   : 2024 district-level 8-month conversion gradient (Inf if no 2024)
  #   sero_err      : population-weighted regional cumulative exposure
  #   sero_dist_err : survey-matched district seroprevalence
  #   rank_pen      : whether 2023 districts are ordered correctly (0/1/2/3)
  #   corr_pen      : ward-level 2023 spatial agreement, (1 - r)^2
  # Names kept short for backward-compat with downstream column refs;
  # semantics are now consistently probability-scale.
  foi_err <- relative_mse(model_regional_seroconv_2023,
                          calibration_targets$obs_regional_seroconv_2023)
  sero_err <- relative_mse(model_sero, calibration_targets$obs_sero)
  sero_sample_err <- relative_mse(model_sero_sample, calibration_targets$obs_sero_sample)
  dist_err <- relative_mse(district_err_df$model_prob, district_err_df$obs)
  sero_dist_err <- relative_mse(district_sero_err_df$model_sero,
                                district_sero_err_df$obs_sero)

  # 2024 components — same estimator (simple mean across surveyed wards) as
  # 2023. If prob_2024_8mo is not supplied, both default to Inf so the joint
  # objective fails loudly rather than silently fitting only 2023.
  has_2024 <- !is.null(prob_2024_8mo_surveyed) && nrow(prob_2024_8mo_surveyed) > 0
  foi_err_24 <- if (has_2024) {
    relative_mse(model_regional_seroconv_2024_8mo,
                 calibration_targets$obs_regional_seroconv_2024_8mo)
  } else {
    Inf
  }
  district_err_df_24 <- calibration_targets$obs_district_seroconv_2024_8mo |>
    left_join(model_district_2024_8mo |> rename(model_prob = prob), by = "district")
  dist_err_24 <- if (has_2024 && nrow(district_err_df_24) > 0) {
    relative_mse(district_err_df_24$model_prob, district_err_df_24$obs)
  } else {
    Inf
  }

  model_rank <- model_district_2023$district
  target_rank <- calibration_targets$obs_district_seroconv_2023$district
  rank_compare <- model_rank[seq_along(target_rank)]
  rank_pen <- sum(is.na(rank_compare) | target_rank != rank_compare)
  corr_pen <- (1 - corr)^2

  # Helper to expose district-specific values in a stable format for saving to
  # calibration output tables and later ABC summaries.
  district_value <- function(df, district_name, value_col) {
    value <- df[[value_col]][df$district == district_name]
    if (length(value) == 0) {
      return(NA_real_)
    }
    value[[1]]
  }

  # ----------------------------------------------------------------------------
  # Objective hierarchy
  # ----------------------------------------------------------------------------
  # Historical objective hierarchy for Paul/Jess probability-scale diagnostics.
  # The final direct-count calibration uses obj_direct_counts, calculated by
  # score_direct_count_run(); these derived objectives are still saved so the
  # final posterior can be compared with the catalytic FOI layer.
  #
  # The five preceding objectives (obj, obj_reg_sero, obj_dist_sero,
  # obj_dist_sero_no_rank, obj_balanced_dist_sero, obj_equal_blocks) are
  # retained because older calibration grids were ranked by them and
  # run_abc_posterior.R's COMPARE mode can dump marginals under each.
  # Computing them here is cheap
  # (single arithmetic on values already in scope) so they cost ~nothing.
  obj_reg_sero <- 6 * foi_err + 4 * dist_err + 2 * sero_err + rank_pen + 2 * corr_pen
  obj_dist_sero <- 6 * foi_err + 4 * dist_err + 2 * sero_dist_err + rank_pen + 2 * corr_pen
  obj_dist_sero_no_rank <- 6 * foi_err + 4 * dist_err + 2 * sero_dist_err + 2 * corr_pen
  obj_balanced_dist_sero <- 3 * foi_err + 3 * dist_err + 2 * sero_dist_err + 2 * corr_pen
  obj_equal_blocks <- foi_err + dist_err + sero_dist_err + corr_pen

  # Within each year the regional:district ratio is 3:2, matching the 6:4
  # ratio of obj_dist_sero_no_rank. No rank penalty because the 2024 rank
  # target is not statistically robust (Babati-Hai ordering is sensitive to
  # one or two high-leverage wards; see 2024 robustness audit).
  obj_joint_2023_2024 <- 3 * foi_err + 2 * dist_err +
                         3 * foi_err_24 + 2 * dist_err_24 +
                         2 * sero_dist_err + 2 * corr_pen

  ward_diagnostics <- c(
    named_metric_values(
      prob_2023_surveyed, "ward", "annual_foi_prob",
      "ward_seroconv_2023_", surveyed_wards
    ),
    named_metric_values(
      prob_2024_8mo_surveyed, "ward", "annual_foi_prob",
      "ward_seroconv_2024_8mo_", surveyed_wards
    ),
    named_metric_values(
      sero_surveyed, "ward_id", "seroprev",
      "ward_sero23_", sero_wards
    )
  )

  list(
    # Default objective: see obj_joint_2023_2024 for the production weight set.
    obj                    = obj_dist_sero,
    obj_dist_sero          = obj_dist_sero,
    obj_reg_sero           = obj_reg_sero,
    obj_dist_sero_no_rank  = obj_dist_sero_no_rank,
    obj_balanced_dist_sero = obj_balanced_dist_sero,
    obj_equal_blocks       = obj_equal_blocks,
    obj_joint_2023_2024    = obj_joint_2023_2024,
    seroconv_2023     = model_regional_seroconv_2023,
    seroconv_2024_8mo = model_regional_seroconv_2024_8mo,
    sero            = model_sero,
    sero_pop        = model_sero,
    sero_sample     = model_sero_sample,
    corr            = corr,
    foi_err         = foi_err,
    foi_err_24      = foi_err_24,
    dist_err        = dist_err,
    dist_err_24     = dist_err_24,
    sero_err        = sero_err,
    sero_reg_err    = sero_err,
    sero_sample_err = sero_sample_err,
    sero_dist_err   = sero_dist_err,
    corr_pen        = corr_pen,
    rank_pen        = rank_pen,
    rank = paste(model_district_2023$district, collapse = ">"),
    babati_seroconv_2023 = district_value(model_district_2023, "babati", "prob"),
    hai_seroconv_2023    = district_value(model_district_2023, "hai",    "prob"),
    moshi_seroconv_2023  = district_value(model_district_2023, "moshi",  "prob"),
    babati_seroconv_2024_8mo = district_value(model_district_2024_8mo, "babati", "prob"),
    hai_seroconv_2024_8mo    = district_value(model_district_2024_8mo, "hai",    "prob"),
    moshi_seroconv_2024_8mo  = district_value(model_district_2024_8mo, "moshi",  "prob"),
    babati_sero23     = district_value(model_district_sero, "babati", "sero"),
    hai_sero23        = district_value(model_district_sero, "hai",    "sero"),
    moshi_sero23      = district_value(model_district_sero, "moshi",  "sero"),
    babati_sero23_pop = district_value(model_district_sero, "babati", "sero_pop"),
    hai_sero23_pop    = district_value(model_district_sero, "hai",    "sero_pop"),
    moshi_sero23_pop  = district_value(model_district_sero, "moshi",  "sero_pop")
  ) |>
    c(ward_diagnostics)
}

#' Extract ward-level seroprevalence at specific model days
#'
#' @param out_df Output from run_ode (data frame with state columns and date)
#' @param ward_ids Ward identifiers [n_wards]
#' @param target_days Numeric vector of model days at which to extract
#' @return Data frame with ward_idx, day, seroprevalence
extract_ward_seroprevalence <- function(out_df, ward_ids, target_days) {
  Sh <- extract_comp_mat(out_df, "Sh", ward_ids)
  Eh <- extract_comp_mat(out_df, "Eh", ward_ids)
  Ih <- extract_comp_mat(out_df, "Ih", ward_ids)
  Rh <- extract_comp_mat(out_df, "Rh", ward_ids)

  results <- vector("list", length(target_days))
  for (k in seq_along(target_days)) {
    row_idx <- which.min(abs(out_df$time - target_days[k]))
    Nh <- Sh[row_idx, ] + Eh[row_idx, ] + Ih[row_idx, ] + Rh[row_idx, ]
    sp <- Rh[row_idx, ] / pmax(Nh, 1)

    results[[k]] <- tibble(
      ward_id   = ward_ids,
      ward_idx  = seq_along(ward_ids),
      model_day = target_days[k],
      seroprev  = as.numeric(sp)
    )
  }
  bind_rows(results)
}
