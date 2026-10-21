# ==============================================================================
# 01_data_prep.R
#
# Data loading and processing: climate, wards, movement, contiguity
#
# Outputs used by the ODE model:
#   - ward_ids, n_wards, N_star, ward_district : ward identifiers and properties
#   - MAT_A_*, MAT_C_*, MAT_TEMP_*            : climate matrices [n_days x n_wards]
#   - MAT_HABITAT_A_*, MAT_HABITAT_C_*        : habitat matrices (currently uniform)
#   - move_mats                                : 12 monthly movement matrices
#   - W_adj, adj_deg                           : adjacency for local coupling
#   - params_static                            : parameter list (VPH wired per district)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(zoo)
  library(Matrix)
  library(igraph)
  library(readr)
})

source("R/00_config.R")

set.seed(422)

# ==============================================================================
# Helper functions
# ==============================================================================

# Standardise ward identifiers to lowercase + trimmed whitespace
clean_id <- function(x) tolower(trimws(x))

# Create named state variable identifiers: e.g. nm("Sh", "babati/ward1") -> "Sh__babati/ward1"
# Used by deSolve to track compartments across wards in a single state vector
nm <- function(prefix, ids) paste0(prefix, "__", ids)

# Recycle a matrix M to fill T_days rows using modular indexing.
# Used to create burn-in climate by cycling observed data multiple times.
# E.g., 25 years of climate data (9131 rows) cycled to fill 75 years (27393 rows)
cycle_to_days <- function(M, T_days) {
  stopifnot(is.matrix(M), nrow(M) > 0)
  idx <- ((seq_len(T_days) - 1) %% nrow(M)) + 1
  M[idx, , drop = FALSE]
}


# ==============================================================================
# Climate-driver matrix builders
# ==============================================================================
#
# These functions translate daily-resolution climate data into the per-ward
# daily recruitment-scaling matrices that the ODE consumes.
#
#   Aedes (standing-water memory): exponential decay of daily rainfall, with
#     district-specific half-life. The water stock carries state across days
#     (and across the burn-in -> experiment boundary) and is mapped to a
#     triangular recruitment curve via scale_aedes_mat().
#
#   Culex (14-day rolling rainfall): fixed-window memory kernel applied via
#     scale_culex_mat(). No internal state.
#
# All three functions accept matrices [n_days x n_wards] so they can be used
# for both the burn-in and experiment periods, and inside per-sample
# calibration loops where the half-life is a fitted parameter.

#' Build a per-ward exponential water-stock from daily rainfall.
#'
#' W_t = R_t + rho * W_{t-1}, with rho = exp(ln(0.5) / h). The half-life h
#' (days) can be a scalar or a per-ward vector. The optional w_max argument
#' caps the stock at a saturation level (e.g. 200 mm-equivalent) to reflect
#' finite breeding-site capacity. The optional initial argument lets the
#' burn-in stock carry into the experiment period.
#'
#' @param rain_mat Daily rainfall matrix [n_days x n_wards]
#' @param half_life_days Scalar or per-ward vector of half-lives (days)
#' @param initial Optional initial-stock vector [n_wards]; defaults to zero
#' @param w_max Optional saturation cap (single value); NULL disables capping
#' @return Matrix [n_days x n_wards] of water-stock values
standing_water_mat <- function(rain_mat, half_life_days,
                               initial = NULL, w_max = NULL) {
  if (length(half_life_days) == 1L) {
    half_life_days <- rep(half_life_days, ncol(rain_mat))
  }
  if (length(half_life_days) != ncol(rain_mat)) {
    stop("half_life_days must be scalar or match rain_mat columns.")
  }
  rho <- exp(log(0.5) / half_life_days)
  prev <- if (is.null(initial)) rep(0, ncol(rain_mat)) else initial
  out <- matrix(0, nrow = nrow(rain_mat), ncol = ncol(rain_mat),
                dimnames = dimnames(rain_mat))
  for (i in seq_len(nrow(rain_mat))) {
    prev <- rain_mat[i, ] + rho * prev
    if (!is.null(w_max)) {
      prev <- pmin(prev, w_max)
    }
    out[i, ] <- prev
  }
  out
}

#' Apply the Aedes triangular recruitment curve to a driver matrix.
#'
#' Maps the water-stock driver (or any analogous matrix) to a scaling factor
#' in [base_frac, 1]. With r_flush > r_peak, the curve ramps up to r_peak and
#' then descends back toward the baseline at r_flush. If the water-stock cap is
#' below r_flush, the stock cannot exceed that cap; the descending limb still
#' applies whenever the capped stock is above r_peak.
#'
#' r_threshold / r_peak / r_flush may each be a scalar or a per-ward vector
#' (length = ncol(driver_mat)). Vector inputs let curve thresholds vary by
#' ward or district.
#'
#' @param driver_mat Driver matrix [n_days x n_wards] (e.g. water stock)
#' @param r_flush Descent-endpoint, scalar or per-ward vector
#' @param r_threshold Ramp-up start, scalar or per-ward vector
#' @param r_peak Ramp-up end / descent start, scalar or per-ward vector
#' @param base_frac_A Dry-season floor (scalar)
#' @return Matrix [n_days x n_wards] of recruitment-scaling values
scale_aedes_mat <- function(driver_mat, r_flush = 200, base_frac_A = 0.05,
                            r_threshold = 50, r_peak = 130) {
  # Broadcast scalar or per-ward vector params to matrices the same shape as
  # driver_mat, so the ifelse() arithmetic is uniformly elementwise.
  expand_to_mat <- function(x, name) {
    if (length(x) == 1L) {
      return(unname(x))
    }
    if (length(x) != ncol(driver_mat)) {
      stop(sprintf("%s length (%d) must be 1 or ncol(driver_mat) (%d)",
                   name, length(x), ncol(driver_mat)))
    }
    # Strip names on x so rep() doesn't drag names into the matrix cells and
    # confuse downstream arithmetic. dimnames stay tied to the driver.
    matrix(rep(unname(x), each = nrow(driver_mat)),
           nrow = nrow(driver_mat), ncol = ncol(driver_mat),
           dimnames = dimnames(driver_mat))
  }
  r_flush_m     <- expand_to_mat(r_flush,     "r_flush")
  r_threshold_m <- expand_to_mat(r_threshold, "r_threshold")
  r_peak_m      <- expand_to_mat(r_peak,      "r_peak")

  if (any(r_flush_m <= r_peak_m)) {
    stop("r_flush must exceed r_peak for the triangular Aedes curve.")
  }
  response <- ifelse(driver_mat < r_threshold_m, 0,
    ifelse(driver_mat <= r_peak_m,
      (driver_mat - r_threshold_m) / (r_peak_m - r_threshold_m),
      pmax(0, (r_flush_m - driver_mat) / (r_flush_m - r_peak_m))
    )
  )
  base_frac_A + (1 - base_frac_A) * response
}

#' Apply the Culex parabolic-window recruitment curve to a rainfall matrix.
#'
#' Wraps calc_culex_recruitment() (from 00_config.R) over a matrix input.
#' Input is 14-day cumulative rainfall (mm).
#'
#' @param rain_mat 14-day cumulative rainfall matrix [n_days x n_wards]
#' @param base_frac_C Dry-season floor (higher than Aedes because Culex breed
#'   in semi-permanent water)
#' @return Matrix [n_days x n_wards] of recruitment-scaling values
scale_culex_mat <- function(rain_mat, base_frac_C = 0.40,
                            w_min = 5, w_max_culex = 120) {
  # Broadcast w_min / w_max_culex (scalar or per-ward vector) to matrices.
  expand_to_mat <- function(x, name) {
    if (length(x) == 1L) {
      return(unname(x))
    }
    if (length(x) != ncol(rain_mat)) {
      stop(sprintf("%s length (%d) must be 1 or ncol(rain_mat) (%d)",
                   name, length(x), ncol(rain_mat)))
    }
    matrix(rep(unname(x), each = nrow(rain_mat)),
           nrow = nrow(rain_mat), ncol = ncol(rain_mat),
           dimnames = dimnames(rain_mat))
  }
  w_min_m   <- expand_to_mat(w_min,       "w_min")
  w_max_m   <- expand_to_mat(w_max_culex, "w_max_culex")

  response <- ifelse(rain_mat < w_min_m | rain_mat > w_max_m, 0,
    4 * (rain_mat - w_min_m) * (w_max_m - rain_mat) / (w_max_m - w_min_m)^2)
  response[response < 0] <- 0
  base_frac_C + (1 - base_frac_C) * response
}

#' End-to-end builder: daily rainfall (burn-in + experiment) -> Aedes
#' recruitment scaling matrices via the standing-water mechanism, with
#' burn-in stock carrying into the experiment.
#'
#' @param rain_daily_burn Daily-rainfall burn-in matrix [T_BURNIN x n_wards]
#' @param rain_daily_exp Daily-rainfall experiment matrix [n_exp_days x n_wards]
#' @param half_lives_per_ward Per-ward half-life vector (days) [n_wards]
#' @param w_max Saturation cap for the standing-water stock
#' @param r_flush Descent endpoint for the Aedes triangular recruitment curve
#' @param base_frac_A Aedes dry-season floor
#' @return List with elements MAT_A_BURN, MAT_A_EXP, stock_burn, stock_exp
build_aedes_water_driver <- function(rain_daily_burn, rain_daily_exp,
                                     half_lives_per_ward,
                                     w_max = 200, r_flush = w_max,
                                     base_frac_A = 0.05) {
  stock_burn <- standing_water_mat(rain_daily_burn, half_lives_per_ward,
                                   w_max = w_max)
  stock_exp <- standing_water_mat(rain_daily_exp, half_lives_per_ward,
                                  initial = stock_burn[nrow(stock_burn), ],
                                  w_max = w_max)
  list(
    MAT_A_BURN = scale_aedes_mat(stock_burn, r_flush = r_flush,
                                 base_frac_A = base_frac_A),
    MAT_A_EXP  = scale_aedes_mat(stock_exp,  r_flush = r_flush,
                                 base_frac_A = base_frac_A),
    stock_burn = stock_burn,
    stock_exp  = stock_exp
  )
}

# ==============================================================================
# Load and process climate data
# ==============================================================================

message("Loading climate data...\n")
daily_climate_df <- readRDS(PATH_CLIMATE)

# Process raw climate into per-ward daily time series with derived columns:
#   rain_3d / rain_14d  : rolling sums over 3 and 14 days
#   scale_aedes/culex   : rainfall-dependent recruitment scaling factors
#                         (see calc_recruit_scale() in 00_config.R)
# Babati Urban wards excluded
ward_climate <- daily_climate_df |>
  mutate(fullname = clean_id(fullname), date = as.Date(date)) |>
  filter(!str_detect(fullname, "babati.*urban")) |>
  group_by(fullname, date) |>
  summarise(
    rain = mean(rain, na.rm = TRUE),
    temp = mean(temp, na.rm = TRUE), .groups = "drop"
  ) |>
  arrange(fullname, date) |>
  group_by(fullname) |>
  mutate(
    rain_3d = zoo::rollsum(rain, k = 3, fill = 0, align = "right"),
    rain_14d = zoo::rollsum(rain, k = 14, fill = 0, align = "right"),
    scale_aedes = calc_recruit_scale(rain_14d, "aedes"),
    scale_culex = calc_recruit_scale(rain_14d, "culex")
  ) |>
  ungroup()

has_temp_data <- "temp" %in% names(daily_climate_df) &&
  sum(!is.na(daily_climate_df$temp)) > 0

message("  Temperature data available:", has_temp_data, "\n")
message(sprintf("  Climate date range: %s to %s",
                min(ward_climate$date), max(ward_climate$date)))

# ==============================================================================
# Load and process ward data
# ==============================================================================

message("\nLoading ward data...\n")
# Ward-level static data: cattle population, elevation, building density
# N_star is the carrying capacity (cattle population) per ward, floored at 1000
# to avoid near-zero denominators in the ODE system
ward_info <- readRDS(PATH_WARD_SUM) |>
  mutate(
    clean_name = clean_id(fullname),
    district = str_split_i(clean_name, "/", 2)
  ) |>
  filter(!str_detect(clean_name, "babati.*urban")) |>
  mutate(N_star = pmax(Cattle_pop, 1000))

# Named vectors for fast lookup by ward ID throughout the model
ward_ids <- sort(ward_info$clean_name)
n_wards <- length(ward_ids)
N_star <- setNames(ward_info$N_star, ward_info$clean_name)[ward_ids]
ward_district <- setNames(ward_info$district, ward_info$clean_name)[ward_ids]

# Build ward-level VPH from district × species values (00_config.R)
# Assign per-ward vectors-per-host (VPH) by district and species.
# The 4 calibration parameters (VPH_A_KILI, VPH_C_KILI, VPH_A_BABATI,
# VPH_C_BABATI) are defined in 00_config.R. Multiplying by mortality rate
# (mu) gives the maximum recruitment rate (lam_max) used in the vector ODE:
#   Lambda = lam_max * N_star * f(rainfall)
vph_a_ward <- ifelse(ward_district == "babati", VPH_A_BABATI, VPH_A_KILI)
vph_c_ward <- ifelse(ward_district == "babati", VPH_C_BABATI, VPH_C_KILI)
names(vph_a_ward) <- names(vph_c_ward) <- ward_ids
params_static$lam_A_max <- params_static$mu_A * vph_a_ward
params_static$lam_C_max <- params_static$mu_C * vph_c_ward

message(sprintf("  Wards: %d  |  Total cattle: %s",
                n_wards, format(sum(N_star), big.mark = ",")))
message("  Districts: ", paste(unique(ward_district), collapse = ", "))
message(sprintf("  VPH: Kili A=%g C=%g | Babati A=%g C=%g",
                VPH_A_KILI, VPH_C_KILI, VPH_A_BABATI, VPH_C_BABATI))

stopifnot(
  length(ward_ids) > 0, length(N_star) == n_wards,
  all(is.finite(N_star)), !anyDuplicated(ward_ids)
)

# ==============================================================================
# Load movement data
# ==============================================================================

move_mats <- NULL
network_metrics <- NULL

if (USE_MOVEMENT) {
  message("Loading movement data...")
  # sim_edge_list is a list of 100 simulated movement realisations (from simulate_movement.R).
  # Each element is an edge list with columns: month, from, to, count.
  # SIM_ID selects which realisation to use (default: 1).
  sim_edge_list <- readRDS(PATH_MOVEMENT)

  # Extract one realisation and keep only edges between study wards
  edge_df <- sim_edge_list[[SIM_ID]] |>
    mutate(from = clean_id(from), to = clean_id(to), month = as.integer(month)) |>
    filter(count > 0, from %in% ward_ids, to %in% ward_ids)

  # Build 12 monthly sparse movement matrices [n_wards x n_wards]
  # Entry (i,j) = number of cattle moved from ward i to ward j in that month
  node_id <- setNames(seq_along(ward_ids), ward_ids)

  move_mats <- lapply(1:12, function(m) {
    dfm <- edge_df |> filter(month == m)
    if (nrow(dfm) == 0) {
      return(Matrix(0, n_wards, n_wards, sparse = TRUE))
    }
    sparseMatrix(
      i = node_id[dfm$from], j = node_id[dfm$to], x = dfm$count,
      dims = c(n_wards, n_wards)
    )
  })

  # Build igraph network for summary metrics (not used by ODE, for diagnostics)
  agg_edges <- edge_df |>
    group_by(from, to) |>
    summarise(weight = sum(count), .groups = "drop")
  g <- graph_from_data_frame(agg_edges, directed = TRUE, vertices = ward_ids)

  network_metrics <- tibble(ward = ward_ids) |>
    mutate(
      out_strength = strength(g, mode = "out"), # total cattle exported
      in_strength = strength(g, mode = "in"), # total cattle imported
      total_strength = strength(g, mode = "all"),
      out_deg = degree(g, mode = "out"), # number of trade partners
      in_deg = degree(g, mode = "in"),
      total_deg = degree(g, mode = "all"),
      betweenness = replace_na(betweenness(g, directed = TRUE), 0)
    )

  message(sprintf("  Total cattle moved (annual): %s | wards with movement: %d",
                  format(sum(edge_df$count), big.mark = ","),
                  sum(network_metrics$total_deg > 0)))
} else {
  message("Market movement disabled.")
  network_metrics <- tibble(
    ward = ward_ids,
    out_deg = 0, in_deg = 0, total_deg = 0, betweenness = 0
  )
}


# ==============================================================================
# Load contiguous ward data (local coupling)
# ==============================================================================

W_adj <- NULL
adj_deg <- NULL

if (USE_LOCAL_COUPLING && file.exists(PATH_CONTIGUOUS)) {
  message("Loading contiguous ward data...")

  contiguous_wards <- read_csv(PATH_CONTIGUOUS, show_col_types = FALSE)

  # Build row-normalised adjacency matrix W_adj for local coupling.
  # W_adj[i,j] = 1/degree(i) if ward j is contiguous with ward i, else 0.
  # This means each ward's contribution from neighbours sums to 1 (diffusion kernel).
  # Symmetrise the edge list since contiguity is undirected.
  cw <- contiguous_wards |>
    mutate(from = clean_id(from), to = clean_id(to)) |>
    filter(from %in% ward_ids, to %in% ward_ids, from != to)
  # Symmetrise: if A borders B, ensure both A->B and B->A exist
  cw <- bind_rows(cw, cw |> transmute(from = to, to = from)) |> distinct(from, to)

  # Map ward names to integer indices for sparse matrix construction
  node_id_adj <- setNames(seq_along(ward_ids), ward_ids)

  # Build binary adjacency matrix A[i,j] = 1 if wards i and j are contiguous
  if (nrow(cw) > 0) {
    A <- sparseMatrix(
      i = node_id_adj[cw$from], j = node_id_adj[cw$to], x = 1,
      dims = c(n_wards, n_wards), dimnames = list(ward_ids, ward_ids)
    )
  } else {
    A <- sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_wards, n_wards), dimnames = list(ward_ids, ward_ids)
    )
  }

  # Row-normalise so each ward's neighbour weights sum to 1:
  # W_adj[i,j] = 1/degree(i) if contiguous, 0 otherwise
  adj_deg <- Matrix::rowSums(A)
  inv_deg <- ifelse(adj_deg > 0, 1 / adj_deg, 0)
  W_adj <- Diagonal(x = inv_deg) %*% A
  dimnames(W_adj) <- list(ward_ids, ward_ids)

  message(sprintf("  Adjacent ward pairs: %d | isolated wards: %d | mean neighbours: %.1f",
                  nrow(cw) / 2, sum(adj_deg == 0),
                  mean(adj_deg[adj_deg > 0])))
} else if (USE_LOCAL_COUPLING) {
  warning("Local coupling enabled but contiguous ward file not found.")
} else {
  message("Local coupling disabled.")
}

# ==============================================================================
# Build climate matrices
# ==============================================================================
#
# The primary matrices are the hindcast period (2000-2024).
# The burn-in cycles these same matrices 3 times — no separate burn-in matrices.

message("Building climate matrices...")

# Pivot a long-format climate data frame into a [n_days x n_wards] matrix.
# Each row is a day, each column is a ward, values are the specified climate variable.
# Missing ward-day combinations are filled with fill_val (defaults: 0.1 for recruitment, 25 for temp).
make_climate_matrix <- function(climate_df, ward_ids, col_name, fill_val = 0.1) {
  climate_df |>
    select(fullname, date, value = all_of(col_name)) |>
    mutate(date = as.Date(date), fullname = clean_id(fullname)) |>
    group_by(fullname, date) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(
      id_cols = date, names_from = fullname,
      values_from = value, values_fill = list(value = fill_val)
    ) |>
    arrange(date) |>
    select(all_of(ward_ids)) |>
    as.matrix()
}

# Keep only study wards for matrix construction
climate_subset <- ward_climate |> filter(fullname %in% ward_ids)

# --- Hindcast matrices: actual 2000-2024 climate (NOT cycled) ---
# Only built when needed (EXPERIMENT_MODE == "hindcast").
# Scripts that only run the contemporary experiment skip this.
need_hindcast <- exists("EXPERIMENT_MODE") && EXPERIMENT_MODE == "hindcast"

if (need_hindcast) {
  climate_hindcast <- climate_subset |>
    filter(date >= HINDCAST_START, date <= HINDCAST_END)

  if (nrow(climate_hindcast) > 0) {
    MAT_A_HINDCAST <- make_climate_matrix(climate_hindcast, ward_ids, "scale_aedes", 0.1)
    MAT_C_HINDCAST <- make_climate_matrix(climate_hindcast, ward_ids, "scale_culex", 0.1)
    MAT_TEMP_HINDCAST <- make_climate_matrix(climate_hindcast, ward_ids, "temp", 25)
    hindcast_dates <- seq.Date(HINDCAST_START, by = "day", length.out = nrow(MAT_A_HINDCAST))

    message(sprintf("  Hindcast: [%d days x %d wards] (%s to %s)",
                    nrow(MAT_A_HINDCAST), n_wards, HINDCAST_START, HINDCAST_END))
  } else {
    stop("No climate data for hindcast period (", HINDCAST_START, " to ", HINDCAST_END, ").")
  }
} else {
  message("  Hindcast matrices: SKIPPED (EXPERIMENT_MODE != 'hindcast')")
}

# --- Contemporary experiment matrices ---
climate_exp <- climate_subset |> filter(date >= EXP_START, date <= EXP_END)

if (nrow(climate_exp) > 0) {
  MAT_A_EXP <- make_climate_matrix(climate_exp, ward_ids, "scale_aedes", 0.1)
  MAT_C_EXP <- make_climate_matrix(climate_exp, ward_ids, "scale_culex", 0.1)
  MAT_TEMP_EXP <- make_climate_matrix(climate_exp, ward_ids, "temp", 25)
  exp_dates <- seq.Date(EXP_START, by = "day", length.out = nrow(MAT_A_EXP))
  message(sprintf("  Contemporary: [%d days x %d wards]", nrow(MAT_A_EXP), n_wards))
} else {
  warning("No climate data for contemporary experiment period. Skipping.")
  MAT_A_EXP <- MAT_C_EXP <- MAT_TEMP_EXP <- NULL
  exp_dates <- NULL
}

# --- Burn-in matrices for contemporary experiment ---
# Filter to the representative climate period for cycling (e.g. 2010-2022).
# This avoids hindcasting during burn-in only
climate_for_cycling <- climate_subset |>
  filter(date >= BURNIN_CLIMATE_START, date <= BURNIN_CLIMATE_END)

if (nrow(climate_for_cycling) == 0) {
  warning("No climate data in burn-in range, falling back to all available data.")
  climate_for_cycling <- climate_subset
}

message(sprintf("  Burn-in climate source: %s to %s",
                BURNIN_CLIMATE_START, BURNIN_CLIMATE_END))

# Build climate matrices from the representative period (e.g. 13 years)
MAT_A_ALL <- make_climate_matrix(climate_for_cycling, ward_ids, "scale_aedes", 0.1)
MAT_C_ALL <- make_climate_matrix(climate_for_cycling, ward_ids, "scale_culex", 0.1)
MAT_TEMP_ALL <- if (has_temp_data) {
  make_climate_matrix(climate_for_cycling, ward_ids, "temp", 25)
} else {
  # No temperature data: use constant 25C
  matrix(25, nrow = nrow(MAT_A_ALL), ncol = n_wards, dimnames = list(NULL, ward_ids))
}

cycle_len <- nrow(MAT_A_ALL)

# Cycle the representative climate to fill the full burn-in duration.
# E.g. 13 years of data repeated ~3 times to cover 40 years.
T_BURNIN <- 365 * EXP_BURNIN_YEARS
MAT_A_BURNIN <- cycle_to_days(MAT_A_ALL, T_BURNIN)
MAT_C_BURNIN <- cycle_to_days(MAT_C_ALL, T_BURNIN)
MAT_TEMP_BURNIN <- cycle_to_days(MAT_TEMP_ALL, T_BURNIN)

# --- Daily rainfall matrices (for the Aedes standing-water driver) ---
# These hold the raw daily rainfall (not pre-scaled), which downstream
# scripts feed into standing_water_mat() / build_aedes_water_driver() with
# posterior-median half-lives to produce per-ward Aedes recruitment scaling.
# Culex still uses the pre-computed scale_culex matrices above (14-day window).
MAT_RAIN_DAILY_ALL <- make_climate_matrix(
  climate_for_cycling, ward_ids, "rain", 0
)
MAT_RAIN_DAILY_BURNIN <- cycle_to_days(MAT_RAIN_DAILY_ALL, T_BURNIN)

if (!is.null(exp_dates)) {
  MAT_RAIN_DAILY_EXP <- make_climate_matrix(climate_exp, ward_ids, "rain", 0)
} else {
  MAT_RAIN_DAILY_EXP <- NULL
}

message(sprintf("  Contemporary burn-in: [%d days x %d wards]", T_BURNIN, n_wards))

# ==============================================================================
# Habitat suitability (uniform — disabled)
# ==============================================================================
# All spatial FOI variation emerges from:
#   - Rainfall differences between wards (recruitment timing)
#   - Temperature differences (vector biology)
#   - Cattle population structure (N_star)
#   - Movement network connectivity
#   - Local coupling (adjacent ward mixing)
# See archive/ for previous INLA and NDVI approaches.

message("  Habitat multipliers: UNIFORM (disabled)")

if (need_hindcast) {
  MAT_HABITAT_A_HINDCAST <- matrix(1,
    nrow = length(hindcast_dates), ncol = n_wards,
    dimnames = list(NULL, ward_ids)
  )
  MAT_HABITAT_C_HINDCAST <- MAT_HABITAT_A_HINDCAST
}

if (!is.null(exp_dates)) {
  MAT_HABITAT_A_EXP <- matrix(1,
    nrow = length(exp_dates), ncol = n_wards,
    dimnames = list(NULL, ward_ids)
  )
  MAT_HABITAT_C_EXP <- MAT_HABITAT_A_EXP
} else {
  MAT_HABITAT_A_EXP <- MAT_HABITAT_C_EXP <- NULL
}

MAT_HABITAT_A_BURNIN <- matrix(1,
  nrow = T_BURNIN, ncol = n_wards,
  dimnames = list(NULL, ward_ids)
)
MAT_HABITAT_C_BURNIN <- MAT_HABITAT_A_BURNIN

# ==============================================================================
# Summary
# ==============================================================================

message("Data preparation complete")
message(sprintf("  Wards: %d", n_wards))
if (need_hindcast) {
  message(sprintf("  Hindcast: %d days (%s to %s)",
                  length(hindcast_dates), HINDCAST_START, HINDCAST_END))
  message(sprintf("  Hindcast burn-in: %d cycles x %d days = %d days",
                  N_BURNIN_CYCLES, length(hindcast_dates),
                  N_BURNIN_CYCLES * length(hindcast_dates)))
}
if (!is.null(exp_dates)) {
  message(sprintf("  Contemporary experiment: %d days", length(exp_dates)))
  message(sprintf("  Contemporary burn-in: %d days (%d years, cycling %s to %s)",
                  T_BURNIN, EXP_BURNIN_YEARS, BURNIN_CLIMATE_START, BURNIN_CLIMATE_END))
}
