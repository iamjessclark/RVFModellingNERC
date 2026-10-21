# ==============================================================================
# 02_model.R
#
# ODE system, movement events, and simulation runners for RVF model
#
# Model structure:
#   Hosts (cattle): SEIRS model with birth/death and waning immunity (constant population N_star)
#   Vectors (Aedes, Culex): SEI model with climate-driven recruitment
#
# State variables per ward (13 total):
#   Sh, Eh, Ih, Rh: Susceptible, Exposed, Infectious, Recovered hosts
#   SvA, EvA, IvA: Susceptible, Exposed, Infectious Aedes vectors
#   SvC, EvC, IvC: Susceptible, Exposed, Infectious Culex vectors
#   Ch:            Cumulative host incidence (for outbreak metrics)
#   Hh_A, Hh_C:    Per-genus cumulative FOI hazard on hosts
#                  (total host hazard Hh = Hh_A + Hh_C; recovered post-hoc when needed)
#
# Coupling mechanisms:
#   1. Market movement: Monthly cattle trade between wards (discrete events)
#   2. Local coupling: Cattle activity-space mixing near ward boundaries (continuous)
#      - Vectors are stationary; cattle move to shared grazing/watering areas
#      - Host->Vector: vectors bite cattle from neighbouring wards
#      - Vector->Host: cattle experience biting hazard from neighbouring wards
#
# Key matrices [dimensions]:
#   MAT_A, MAT_C: Recruitment scaling [n_days x n_wards] - temporal (rainfall)
#   MAT_HABITAT_A, MAT_HABITAT_C: Habitat suitability [n_days x n_wards] - spatial (currently uniform)
#   MAT_TEMP: Temperature [n_days x n_wards] - for temperature-dependent parameters
#   W_adj: Row-normalised adjacency [n_wards x n_wards] - for local coupling
#   move_mats: List of 12 monthly movement matrices [n_wards x n_wards each]
# ==============================================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(Matrix)
  library(tibble)
})

# Local fallback for state-name construction.
# This helper is usually defined in 01_data_prep.R, but several functions in
# this file depend on it directly. Defining it here as well makes 02_model.R
# safe to source or use independently in interactive sessions.
if (!exists("nm", mode = "function")) {
  nm <- function(prefix, ids) paste0(prefix, "__", ids)
}

# ==============================================================================
# ODE system
# ==============================================================================
# The ODE function is called at each time step by deSolve.
# It receives the current state vector y and returns derivatives dy/dt.
#
# Recruitment formula:
#   Lambda = lam_max * N_star * scale * habitat
#   where:
#     lam_max = base recruitment rate (from VPH_MAX * mu)
#     N_star  = target cattle population (vectors scale with hosts)
#     scale   = temporal variation from rainfall (0-1)
#     habitat = spatial multiplier (currently uniform = 1 for all wards)
#
# Local coupling interpretation:
#   eps = fraction of time cattle spend near ward boundaries
#   Vectors don't move; cattle do (shared grazing, watering points)

#' Convert annualised exogenous infection pressure to a daily value
#'
#' intro_rate can be:
#'   - scalar: same value for all wards
#'   - named vector by ward ID
#' Converted to a daily quantity as 1 - (1 - annual_rate)^(1/365).
#' The returned value is added to effective host infectious prevalence seen by
#' vectors; it is not an explicit host movement event.
build_intro_daily <- function(ward_ids, params) {
  intro <- params$intro_rate %||% 0
  if (length(intro) == 1) {
    daily <- 1 - (1 - intro)^(1 / 365)
    return(rep(daily, length(ward_ids)))
  }
  # Named vector by ward ID
  daily_vec <- rep(0, length(ward_ids))
  names(daily_vec) <- ward_ids
  for (w in ward_ids) {
    val <- intro[w]
    if (is.na(val)) val <- 0
    daily_vec[w] <- 1 - (1 - val)^(1 / 365)
  }
  daily_vec
}


#' Extract final state vector from ODE output
#'
#' @param out_df Output from run_ode
#' @return Named state vector from the final time point
extract_final_state <- function(out_df) {
  out_df |>
    select(-date, -time) |>
    slice_tail(n = 1) |>
    unlist(use.names = TRUE)
}


#' Summarise state at the end of a burn-in cycle
#' quantities that measures for endemic persistence and transmission intensity.
summarise_cycle_state <- function(state, ward_ids) {
  tibble(ward = ward_ids) |>
    mutate(
      Sh = state[nm("Sh", ward)],
      Eh = state[nm("Eh", ward)],
      Ih = state[nm("Ih", ward)],
      Rh = state[nm("Rh", ward)],
      SvA = state[nm("SvA", ward)],
      EvA = state[nm("EvA", ward)],
      IvA = state[nm("IvA", ward)],
      SvC = state[nm("SvC", ward)],
      EvC = state[nm("EvC", ward)],
      IvC = state[nm("IvC", ward)],
      Nh = pmax(Sh + Eh + Ih + Rh, 1),
      prevalence = Ih / Nh,
      seroprevalence = Rh / Nh,
      aedes_vph = (SvA + EvA + IvA) / Nh,
      culex_vph = (SvC + EvC + IvC) / Nh
    ) |>
    select(ward, prevalence, seroprevalence, aedes_vph, culex_vph)
}


#' Compare successive cycle-end states
#'
#' @param prev_state Cycle-end state from the previous cycle
#' @param curr_state Cycle-end state from the current cycle
#' @param ward_ids Ward identifiers
#' @return Scalar convergence metric (smaller = closer to periodic equilibrium)
calc_cycle_delta <- function(prev_state, curr_state, ward_ids) {
  prev_summary <- summarise_cycle_state(prev_state, ward_ids)
  curr_summary <- summarise_cycle_state(curr_state, ward_ids)

  prev_summary |>
    left_join(curr_summary, by = "ward", suffix = c("_prev", "_curr")) |>
    summarise(
      delta_prev = max(abs(prevalence_curr - prevalence_prev), na.rm = TRUE),
      delta_sero = max(abs(seroprevalence_curr - seroprevalence_prev), na.rm = TRUE),
      delta_aedes = max(abs(aedes_vph_curr - aedes_vph_prev), na.rm = TRUE),
      delta_culex = max(abs(culex_vph_curr - culex_vph_prev), na.rm = TRUE)
    ) |>
    unlist(use.names = FALSE) |>
    max(na.rm = TRUE)
}


#' Run burn-in until the seasonal cycle stabilises
#' The system is seasonally forced, so the target is a repeating annual or
#' multi-year cycle rather than a flat equilibrium. We therefore compare the
#' state at the same point in successive climate cycles.
run_burnin_to_periodic <- function(MAT_A, MAT_C, dates, y_init, params,
                                   MAT_TEMP = NULL,
                                   use_movement = TRUE, move_mats = NULL,
                                   move_freq = c("monthly", "fixed_step"),
                                   move_step_days = 14,
                                   move_scale_to_time = TRUE,
                                   move_mult = 1,
                                   ward_ids, N_star,
                                   rebalance_hosts = TRUE,
                                   rebalance_method = "scale",
                                   use_temp_dependence = FALSE,
                                   use_local_coupling = FALSE,
                                   W_adj = NULL, adj_deg = NULL,
                                   MAT_HABITAT_A = NULL, MAT_HABITAT_C = NULL,
                                   min_cycles = 2,
                                   max_cycles = 6,
                                   tol = 0.01,
                                   consecutive_hits = 1,
                                   verbose = FALSE) {
  move_freq <- match.arg(move_freq)

  current_state <- y_init
  prev_end_state <- NULL
  stable_hits <- 0
  cycle_history <- vector("list", max_cycles)
  last_out <- NULL

  for (cycle_id in seq_len(max_cycles)) {
    out_cycle <- run_ode(
      MAT_A = MAT_A,
      MAT_C = MAT_C,
      MAT_TEMP = MAT_TEMP,
      MAT_HABITAT_A = MAT_HABITAT_A,
      MAT_HABITAT_C = MAT_HABITAT_C,
      dates = dates,
      y_init = current_state,
      params = params,
      use_movement = use_movement,
      move_mats = move_mats,
      move_freq = move_freq,
      move_step_days = move_step_days,
      move_scale_to_time = move_scale_to_time,
      move_mult = move_mult,
      ward_ids = ward_ids,
      N_star = N_star,
      rebalance_hosts = rebalance_hosts,
      rebalance_method = rebalance_method,
      use_temp_dependence = use_temp_dependence,
      use_local_coupling = use_local_coupling,
      W_adj = W_adj,
      adj_deg = adj_deg,
      MAT_HABITAT_A = MAT_HABITAT_A,
      MAT_HABITAT_C = MAT_HABITAT_C
    )

    end_state <- extract_final_state(out_cycle) |>
      reset_cumulative(ward_ids = ward_ids)

    cycle_delta <- if (is.null(prev_end_state)) {
      NA_real_
    } else {
      calc_cycle_delta(prev_end_state, end_state, ward_ids)
    }

    cycle_history[[cycle_id]] <- tibble(
      cycle = cycle_id,
      delta = cycle_delta
    )

    if (!is.null(prev_end_state) &&
      cycle_id >= min_cycles &&
      is.finite(cycle_delta) &&
      cycle_delta < tol) {
      stable_hits <- stable_hits + 1
    } else {
      stable_hits <- 0
    }

    if (verbose) {
      cat(
        "    Burn-in cycle", cycle_id,
        "| delta =", ifelse(is.na(cycle_delta), "NA", signif(cycle_delta, 3)),
        "| stable hits =", stable_hits, "\n"
      )
    }

    last_out <- out_cycle
    current_state <- end_state
    prev_end_state <- end_state

    if (stable_hits >= consecutive_hits) {
      break
    }
  }

  list(
    state = current_state,
    out_last = last_out,
    history = bind_rows(cycle_history),
    converged = stable_hits >= consecutive_hits,
    n_cycles = cycle_id
  )
}

#' Create ODE function with all features
#'
#' @param MAT_A Aedes recruitment scaling [n_days x n_wards]
#' @param MAT_C Culex recruitment scaling [n_days x n_wards]
#' @param MAT_TEMP Temperature matrix [n_days x n_wards]
#' @param dates Date sequence [n_days]
#' @param MAT_HABITAT_A Aedes habitat suitability [n_days x n_wards], NULL = uniform
#' @param MAT_HABITAT_C Culex habitat suitability [n_days x n_wards], NULL = uniform
#' @param n_days Number of simulation days
#' @param ward_ids Vector of ward identifiers [n_wards]
#' @param N_star Named vector of cattle populations [n_wards]
#' @param params Parameter list from 00_config.R
#' @param use_temp_dependence Use temperature-dependent vector parameters?
#' @param W_adj Row-normalised adjacency matrix [n_wards x n_wards], NULL = no coupling
#' @param adj_deg Degree vector [n_wards] for handling isolated wards
#' @return ODE function for use with deSolve::ode()
make_ode_coupled <- function(MAT_A, MAT_C, MAT_TEMP, dates, n_days, ward_ids, N_star,
                             params, use_temp_dependence = TRUE,
                             W_adj = NULL, adj_deg = NULL,
                             MAT_HABITAT_A = NULL, MAT_HABITAT_C = NULL) {
  # Pre-compute (evaluated once, not at every time step)
  use_coupling <- !is.null(W_adj)
  eps <- if (use_coupling) params$eps_local else 0
  use_habitat <- !is.null(MAT_HABITAT_A) && !is.null(MAT_HABITAT_C)
  omega_h <- params$omega_h %||% 0 # Waning immunity rate (0 = permanent immunity)
  intro_daily <- build_intro_daily(ward_ids, params)


  # Return the ODE function
  function(t, y, parms) {
    # Convert continuous time to integer day (matrix row index)
    t_int <- max(1, min(n_days, floor(t)))

    # ==========================================================================
    # Look up time-varying inputs from matrices
    # Each lookup returns a vector of length n_wards
    # ==========================================================================

    # Rainfall-driven recruitment scaling (0 to 1)
    scale_A <- MAT_A[t_int, ]
    scale_C <- MAT_C[t_int, ]

    # Habitat suitability multipliers (mean ~1)
    if (use_habitat) {
      habitat_A <- MAT_HABITAT_A[t_int, ]
      habitat_C <- MAT_HABITAT_C[t_int, ]
    } else {
      habitat_A <- rep(1.0, length(ward_ids))
      habitat_C <- rep(1.0, length(ward_ids))
    }

    # Temperature (degrees Celsius)
    temp_vec <- MAT_TEMP[t_int, ]

    # ==========================================================================
    # Extract state variables from y
    # y is a named vector; nm() creates names like "Sh__ward1"
    # ==========================================================================

    # Host compartments [n_wards each]
    Sh <- y[nm("Sh", ward_ids)]
    Eh <- y[nm("Eh", ward_ids)]
    Ih <- y[nm("Ih", ward_ids)]
    Rh <- y[nm("Rh", ward_ids)]
    Nh <- pmax(Sh + Eh + Ih + Rh, 1) # Ensure positive to avoid division by zero

    # Vector compartments [n_wards each]
    SvA <- y[nm("SvA", ward_ids)]
    EvA <- y[nm("EvA", ward_ids)]
    IvA <- y[nm("IvA", ward_ids)]

    SvC <- y[nm("SvC", ward_ids)]
    EvC <- y[nm("EvC", ward_ids)]
    IvC <- y[nm("IvC", ward_ids)]

    # Cumulative trackers [n_wards each] -- not read inside the ODE, only
    # written. Listed here for reference:
    # Ch    : cumulative incidence
    # Hh_A  : cumulative Aedes-mediated host hazard
    # Hh_C  : cumulative Culex-mediated host hazard

    # ==========================================================================
    # Vector parameters (temperature-dependent or static)
    # ==========================================================================

    if (use_temp_dependence) {
      # Use thermal functions from 00_config.R
      # %||% is the null-coalescing operator: returns right side if left is NULL
      # Example: params$temp_species_A %||% "aedes" means
      #          "use params$temp_species_A if not NULL, else use 'aedes'"
      spA <- params$temp_species_A %||% "aedes"
      spC <- params$temp_species_C %||% "culex"

      # Aedes parameters (vary by ward based on local temperature)
      # Minimum floors ensure numerical stability:
      #   - Biting rate: 0.001 (near-zero; mosquitoes don't bite when cold)
      #   - Mortality: 0.01 = 100-day lifespan max (must eventually die)
      #   - EIP rate: 0.01 = 100-day EIP max (virus must eventually develop)
      a_A <- pmax(calc_biting_rate(temp_vec, spA), 0.001)
      mu_A <- pmax(calc_mortality_rate(temp_vec, spA), 0.01)
      sigma_A <- pmax(calc_eip_rate(temp_vec, spA), 0.01)

      # Culex parameters
      a_C <- pmax(calc_biting_rate(temp_vec, spC), 0.001)
      mu_C <- pmax(calc_mortality_rate(temp_vec, spC), 0.01)
      sigma_C <- pmax(calc_eip_rate(temp_vec, spC), 0.01)

      # Transmission probabilities are constant (not temperature-dependent)
      beta_vh_A <- params$beta_vh_A
      beta_vh_C <- params$beta_vh_C
      beta_hv_A <- params$beta_hv_A
      beta_hv_C <- params$beta_hv_C
    } else {
      # Use static parameters (same for all wards and times)
      a_A <- params$a_A
      mu_A <- params$mu_A
      sigma_A <- params$sigma_A
      a_C <- params$a_C
      mu_C <- params$mu_C
      sigma_C <- params$sigma_C
      beta_vh_A <- params$beta_vh_A
      beta_vh_C <- params$beta_vh_C
      beta_hv_A <- params$beta_hv_A
      beta_hv_C <- params$beta_hv_C
    }

    # ==========================================================================
    # Vector recruitment
    # Lambda = lam_max * N_star * scale * habitat
    # Units: vectors/day entering each ward
    # ==========================================================================

    Lambda_A <- params$lam_A_max * N_star * scale_A * habitat_A
    Lambda_C <- params$lam_C_max * N_star * scale_C * habitat_C

    # ==========================================================================
    # Local coupling (bidirectional)
    # Cattle activity-space mixing near ward boundaries
    # Vectors are stationary; cattle move to shared grazing/watering areas
    # eps = fraction of time cattle spend near boundaries (e.g., 0.01 = 1%)
    #
    # Host -> Vector: Vectors in ward i bite cattle that may come from ward j
    #   -> Mix host prevalence (Ih/Nh)
    #
    # Vector -> Host: Cattle from ward i spend time in ward j and get bitten there
    #   -> Mix biting hazard (FOI), not infected vector counts
    # ==========================================================================

    # --- Host -> Vector coupling: mix host prevalence ---
    p_local <- Ih / Nh # Local host prevalence [n_wards]
    p_local_total <- pmin(p_local + intro_daily, 0.95)

    if (use_coupling && eps > 0) {
      # Neighbour-averaged prevalence (W_adj is row-normalised, so this is a weighted mean)
      p_nei <- as.numeric(W_adj %*% p_local_total)

      # Handle isolated wards (no neighbours): use local prevalence
      if (!is.null(adj_deg)) {
        p_nei <- ifelse(adj_deg > 0, p_nei, p_local_total)
      }

      # Coupled prevalence: mix of local and neighbour
      p_coupled <- (1 - eps) * p_local_total + eps * p_nei
    } else {
      p_coupled <- p_local_total
    }

    # --- Vector -> Host coupling: mix biting hazard (FOI), tracked by genus ---
    # Local hazard from each genus separately. Each ward's hazard uses its own
    # Nh denominator (we mix exposure risk between wards, not infected vector
    # counts, because vectors are stationary).
    foi_h_A_local <- a_A * beta_vh_A * IvA / Nh
    foi_h_C_local <- a_C * beta_vh_C * IvC / Nh

    if (use_coupling && eps > 0) {
      # Apply ward-neighbour mixing separately to each genus so that
      # genus-specific cumulative hazards remain meaningful and sum to the
      # total coupled hazard.
      foi_h_A_nei <- as.numeric(W_adj %*% foi_h_A_local)
      foi_h_C_nei <- as.numeric(W_adj %*% foi_h_C_local)

      if (!is.null(adj_deg)) {
        foi_h_A_nei <- ifelse(adj_deg > 0, foi_h_A_nei, foi_h_A_local)
        foi_h_C_nei <- ifelse(adj_deg > 0, foi_h_C_nei, foi_h_C_local)
      }

      foi_h_A <- (1 - eps) * foi_h_A_local + eps * foi_h_A_nei
      foi_h_C <- (1 - eps) * foi_h_C_local + eps * foi_h_C_nei
    } else {
      foi_h_A <- foi_h_A_local
      foi_h_C <- foi_h_C_local
    }

    # Total host hazard (used for incidence). By construction
    # foi_h = foi_h_A + foi_h_C identically, both before and after coupling.
    foi_h <- foi_h_A + foi_h_C

    # --- Force of infection on vectors ---
    # p_coupled includes the low-level exogenous eta_i term, so vectors respond
    # to both within-system host infection and background outside pressure.
    foi_vA <- a_A * beta_hv_A * p_coupled
    foi_vC <- a_C * beta_hv_C * p_coupled

    # ==========================================================================
    # Host dynamics (SEIR with constant population)
    # Birth rate = death rate * N_star maintains constant population
    # ==========================================================================

    B <- params$mu_h * N_star # Birth rate [n_wards]
    inc_h <- foi_h * Sh # New infections [n_wards]

    # ==========================================================================
    # Return derivatives
    # Each setNames() creates named derivatives for deSolve to track
    # ==========================================================================

    list(c(
      # Host derivatives (SEIRS with waning immunity) [n_wards each]
      # dS/dt = births + waning - infection - death
      setNames(B - inc_h - params$mu_h * Sh + omega_h * Rh, nm("Sh", ward_ids)),
      # dE/dt = infection - progression - death
      setNames(inc_h - params$sigma_h * Eh - params$mu_h * Eh, nm("Eh", ward_ids)),
      # dI/dt = progression - recovery - death
      setNames(params$sigma_h * Eh - params$gamma_h * Ih - params$mu_h * Ih, nm("Ih", ward_ids)),
      # dR/dt = recovery - death - waning (omega_h: R->S, 0 = permanent immunity)
      setNames(params$gamma_h * Ih - params$mu_h * Rh - omega_h * Rh, nm("Rh", ward_ids)),

      # Aedes derivatives (SEI) [n_wards each]
      setNames(Lambda_A - foi_vA * SvA - mu_A * SvA, nm("SvA", ward_ids)),
      setNames(foi_vA * SvA - sigma_A * EvA - mu_A * EvA, nm("EvA", ward_ids)),
      setNames(sigma_A * EvA - mu_A * IvA, nm("IvA", ward_ids)),

      # Culex derivatives (SEI) [n_wards each]
      setNames(Lambda_C - foi_vC * SvC - mu_C * SvC, nm("SvC", ward_ids)),
      setNames(foi_vC * SvC - sigma_C * EvC - mu_C * EvC, nm("EvC", ward_ids)),
      setNames(sigma_C * EvC - mu_C * IvC, nm("IvC", ward_ids)),

      # Cumulative trackers [n_wards each]
      # Total cumulative hazard Hh = Hh_A + Hh_C (post-processing if needed).
      setNames(inc_h,   nm("Ch",   ward_ids)),  # cumulative incidence
      setNames(foi_h_A, nm("Hh_A", ward_ids)),  # Aedes-mediated cumulative hazard
      setNames(foi_h_C, nm("Hh_C", ward_ids))   # Culex-mediated cumulative hazard
    ))
  }
}


#' Create simple ODE function (no coupling, no temperature dependence, no habitat)
#' Kept for comparison, debugging, and faster runs
make_ode_static <- function(MAT_A, MAT_C, dates, n_days, ward_ids, N_star, params) {
  intro_daily <- build_intro_daily(ward_ids, params)
  omega_h <- params$omega_h %||% 0


  function(t, y, parms) {
    t_int <- max(1, min(n_days, floor(t)))

    scale_A <- MAT_A[t_int, ]
    scale_C <- MAT_C[t_int, ]

    Sh <- y[nm("Sh", ward_ids)]
    Eh <- y[nm("Eh", ward_ids)]
    Ih <- y[nm("Ih", ward_ids)]
    Rh <- y[nm("Rh", ward_ids)]
    Nh <- pmax(Sh + Eh + Ih + Rh, 1)

    SvA <- y[nm("SvA", ward_ids)]
    EvA <- y[nm("EvA", ward_ids)]
    IvA <- y[nm("IvA", ward_ids)]
    SvC <- y[nm("SvC", ward_ids)]
    EvC <- y[nm("EvC", ward_ids)]
    IvC <- y[nm("IvC", ward_ids)]

    Lambda_A <- params$lam_A_max * N_star * scale_A
    Lambda_C <- params$lam_C_max * N_star * scale_C

    p_local_total <- pmin(Ih / Nh + intro_daily, 0.95)

    # Genus-specific host hazard (no coupling in this static variant).
    foi_h_A <- params$a_A * params$beta_vh_A * IvA / Nh
    foi_h_C <- params$a_C * params$beta_vh_C * IvC / Nh
    foi_h   <- foi_h_A + foi_h_C
    foi_vA <- params$a_A * params$beta_hv_A * p_local_total
    foi_vC <- params$a_C * params$beta_hv_C * p_local_total

    B <- params$mu_h * N_star
    inc_h <- foi_h * Sh

    list(c(
      # dS/dt = births + waning - infection - death
      setNames(B - inc_h - params$mu_h * Sh + omega_h * Rh, nm("Sh", ward_ids)),
      # dE/dt = infection - progression - death
      setNames(inc_h - params$sigma_h * Eh - params$mu_h * Eh, nm("Eh", ward_ids)),
      # dI/dt = progression - recovery - death
      setNames(params$sigma_h * Eh - params$gamma_h * Ih - params$mu_h * Ih, nm("Ih", ward_ids)),
      # dR/dt = recovery - death - waning
      setNames(params$gamma_h * Ih - params$mu_h * Rh - omega_h * Rh, nm("Rh", ward_ids)),
      setNames(Lambda_A - foi_vA * SvA - params$mu_A * SvA, nm("SvA", ward_ids)),
      setNames(foi_vA * SvA - params$sigma_A * EvA - params$mu_A * EvA, nm("EvA", ward_ids)),
      setNames(params$sigma_A * EvA - params$mu_A * IvA, nm("IvA", ward_ids)),
      setNames(Lambda_C - foi_vC * SvC - params$mu_C * SvC, nm("SvC", ward_ids)),
      setNames(foi_vC * SvC - params$sigma_C * EvC - params$mu_C * EvC, nm("EvC", ward_ids)),
      setNames(params$sigma_C * EvC - params$mu_C * IvC, nm("IvC", ward_ids)),
      setNames(inc_h,   nm("Ch",   ward_ids)),
      setNames(foi_h_A, nm("Hh_A", ward_ids)),  # Aedes-mediated cumulative hazard
      setNames(foi_h_C, nm("Hh_C", ward_ids))   # Culex-mediated cumulative hazard
    ))
  }
}


# ==============================================================================
# Movement event functions (market/long-distance trade)
# ==============================================================================
# Market movement is implemented as discrete events (not continuous ODE).
# At each event time, cattle are moved between wards according to the
# movement matrix for that month.
#
# Movement assumptions:
#   - Cattle move proportionally by compartment (infected cattle spread disease)
#   - Outflow is capped at current population (can't move more than you have)
#   - After movement, populations are "rebalanced" to N_star
#
# Rebalancing methods:
#   "scale": Scale all compartments proportionally (preserves S:E:I:R ratios)
#   "susceptible": Add/remove from S only (dilutes immunity when adding)

#' Apply proportional SEIR movement between wards
#'
#' @param S,E,I,R Current compartment sizes [n_wards]
#' @param M Movement matrix [n_wards x n_wards], M[i,j] = cattle from ward i to ward j
#' @param N_star Target populations [n_wards]
#' @param rebalance Apply rebalancing after movement?
#' @param method Rebalancing method: "scale" (recommended) or "susceptible"
#' @return List with updated S, E, I, R vectors [n_wards each]
apply_movement <- function(S, E, I, R, M, N_star, rebalance = TRUE,
                           method = c("scale", "susceptible")) {
  method <- match.arg(method)

  N <- S + E + I + R # Current populations [n_wards]
  O <- Matrix::rowSums(M) # Total outflow from each ward [n_wards]

  # Scale movement if outflow exceeds population (can't move more than you have)
  scale <- ifelse(O > 0, pmin(1, N / O), 0)
  M2 <- Matrix::Diagonal(x = scale) %*% M
  O2 <- Matrix::rowSums(M2)

  # Move each compartment proportionally
  move_comp <- function(X) {
    w <- ifelse(N > 0, X / pmax(N, 1e-9), 0) # Fraction in this compartment
    pmax(X - w * O2 + as.numeric(Matrix::t(M2) %*% w), 0)
  }

  S2 <- move_comp(S)
  E2 <- move_comp(E)
  I2 <- move_comp(I)
  R2 <- move_comp(R)

  # Rebalance to N_star (correct for numerical drift)
  if (rebalance) {
    N2 <- S2 + E2 + I2 + R2

    if (method == "scale") {
      # Scale all compartments proportionally (preserves disease ratios)
      # Treats discrepancy as census correction
      scale_factor <- ifelse(N2 > 0, N_star / N2, 1)
      S2 <- S2 * scale_factor
      E2 <- E2 * scale_factor
      I2 <- I2 * scale_factor
      R2 <- R2 * scale_factor
    } else {
      # All adjustment into S (dilutes immunity when adding cattle)
      delta <- N_star - N2
      S2 <- pmax(S2 + delta, 0)
    }
  }

  list(S = S2, E = E2, I = I2, R = R2)
}


#' Create movement event function for deSolve
#'
#' @param dates Date sequence for simulation [n_days]
#' @param ward_ids Ward identifiers [n_wards]
#' @param move_mats List of 12 monthly movement matrices [n_wards x n_wards each]
#' @param N_star Target populations [n_wards]
#' @param rebalance Apply rebalancing after movement?
#' @param method Rebalancing method
#' @param freq Event frequency: "monthly" (1st of month) or "fixed_step"
#' @param step_days Days between events if freq = "fixed_step"
#' @param scale_to_time Scale movement to time step? (keeps annual totals constant)
#' @param move_mult Movement multiplier for sensitivity analysis
#' @return Event function for deSolve
make_event <- function(dates, ward_ids, move_mats, N_star,
                       rebalance = TRUE,
                       method = c("scale", "susceptible"),
                       freq = c("monthly", "fixed_step"),
                       step_days = 14,
                       scale_to_time = TRUE,
                       move_mult = 1) {
  method <- match.arg(method)
  freq <- match.arg(freq)

  n_days <- length(dates)

  # Determine event times
  if (freq == "monthly") {
    event_times <- which(lubridate::day(dates) == 1)
  } else {
    event_times <- seq(1, n_days, by = step_days)
  }

  function(t, y, parms) {
    ti <- as.integer(round(t))
    if (!(ti %in% event_times)) {
      return(y)
    }

    this_date <- dates[ti]
    m <- lubridate::month(this_date)
    M <- move_mats[[m]]

    # Scale movement to time step if using fixed_step (keeps annual totals constant)
    if (freq == "fixed_step" && isTRUE(scale_to_time)) {
      days_in_month <- lubridate::days_in_month(this_date)
      frac <- step_days / as.numeric(days_in_month)
      M <- M * frac
    }

    # Apply movement multiplier (for sensitivity analysis)
    if (!is.null(move_mult) && is.finite(move_mult) && move_mult != 1) {
      M <- M * move_mult
    }

    # Apply movement to host compartments
    moved <- apply_movement(
      y[nm("Sh", ward_ids)], y[nm("Eh", ward_ids)],
      y[nm("Ih", ward_ids)], y[nm("Rh", ward_ids)],
      M, N_star,
      rebalance = rebalance, method = method
    )

    y[nm("Sh", ward_ids)] <- moved$S
    y[nm("Eh", ward_ids)] <- moved$E
    y[nm("Ih", ward_ids)] <- moved$I
    y[nm("Rh", ward_ids)] <- moved$R
    y
  }
}


# ==============================================================================
# Simulation runner
# ==============================================================================

#' Reset cumulative trackers to zero
#' Call before starting a new experiment from an endemic state
#'
#' @param state State vector from end of burn-in
#' @param ward_ids Ward identifiers
#' @return State vector with Ch, Hh_A, and Hh_C reset to zero
reset_cumulative <- function(state, ward_ids) {
  out <- state
  out[nm("Ch",   ward_ids)] <- 0
  out[nm("Hh_A", ward_ids)] <- 0
  out[nm("Hh_C", ward_ids)] <- 0
  out
}


#' Run ODE simulation with all features
#'
#' Main entry point for running simulations. Sets up the appropriate ODE
#' function and events based on requested features.
#'
#' @param MAT_A Aedes recruitment scaling [n_days x n_wards]
#' @param MAT_C Culex recruitment scaling [n_days x n_wards]
#' @param MAT_TEMP Temperature matrix [n_days x n_wards], NULL = constant 25°C
#' @param MAT_HABITAT_A Aedes habitat suitability [n_days x n_wards], NULL = uniform (1)
#' @param MAT_HABITAT_C Culex habitat suitability [n_days x n_wards], NULL = uniform (1)
#' @param dates Date sequence [n_days]
#' @param y_init Initial state vector [n_wards * 13 compartments]
#' @param params Parameter list from 00_config.R
#' @param use_movement Enable market movement events?
#' @param move_mats List of 12 monthly movement matrices
#' @param move_freq Movement frequency: "monthly" or "fixed_step"
#' @param move_step_days Days between events if move_freq = "fixed_step"
#' @param move_scale_to_time Scale movement to time step?
#' @param move_mult Movement multiplier for sensitivity analysis
#' @param ward_ids Ward identifiers [n_wards]
#' @param N_star Target cattle populations [n_wards]
#' @param rebalance_hosts Rebalance to N_star after movement?
#' @param rebalance_method "scale" (preserve ratios) or "susceptible" (add to S)
#' @param use_temp_dependence Use temperature-dependent vector parameters?
#' @param use_local_coupling Use adjacent ward coupling?
#' @param W_adj Row-normalised adjacency matrix [n_wards x n_wards]
#' @param adj_deg Degree vector [n_wards]
#' @return Data frame with time series of all state variables plus date column
run_ode <- function(MAT_A, MAT_C, dates, y_init, params,
                    MAT_TEMP = NULL,
                    use_movement = TRUE, move_mats = NULL,
                    move_freq = c("monthly", "fixed_step"),
                    move_step_days = 14,
                    move_scale_to_time = TRUE,
                    move_mult = 1,
                    ward_ids, N_star,
                    rebalance_hosts = TRUE,
                    rebalance_method = "scale",
                    use_temp_dependence = FALSE,
                    use_local_coupling = FALSE,
                    W_adj = NULL, adj_deg = NULL,
                    MAT_HABITAT_A = NULL, MAT_HABITAT_C = NULL,
                    ode_method = c("lsoda", "rk4"),
                    rk4_dt = 1) {
  move_freq <- match.arg(move_freq)
  ode_method <- match.arg(ode_method)
  n_days <- length(dates)

  if (ode_method == "rk4") {
    # For fixed-step RK4 we build an explicit sub-daily time grid. Integer days
    # are always included so movement events still occur on the intended dates,
    # and we can return daily summaries in the same shape as the lsoda output.
    times <- sort(unique(round(c(seq(1, n_days, by = rk4_dt), seq_len(n_days)), 10)))
  } else {
    times <- 1:n_days
  }

  # Validate flags

  use_movement <- isTRUE(use_movement) && !is.null(move_mats)
  use_local_coupling <- isTRUE(use_local_coupling) && !is.null(W_adj)

  # Choose ODE function based on features requested
  if (use_temp_dependence || use_local_coupling ||
    !is.null(MAT_HABITAT_A) || !is.null(MAT_HABITAT_C)) {
    # Full-featured ODE (temperature, coupling, and/or habitat)
    if (is.null(MAT_TEMP)) {
      # Default: constant 25°C
      MAT_TEMP <- matrix(25, nrow = n_days, ncol = length(ward_ids))
      colnames(MAT_TEMP) <- ward_ids
    }

    ode_fun <- make_ode_coupled(
      MAT_A, MAT_C, MAT_TEMP, dates, n_days, ward_ids, N_star, params,
      use_temp_dependence = use_temp_dependence,
      W_adj = if (use_local_coupling) W_adj else NULL,
      adj_deg = if (use_local_coupling) adj_deg else NULL,
      MAT_HABITAT_A = MAT_HABITAT_A,
      MAT_HABITAT_C = MAT_HABITAT_C
    )
  } else {
    # Simple ODE (faster, fewer features)
    ode_fun <- make_ode_static(MAT_A, MAT_C, dates, n_days, ward_ids, N_star, params)
  }

  # Set up movement events
  if (use_movement) {
    ev_fun <- make_event(
      dates, ward_ids, move_mats, N_star,
      rebalance = rebalance_hosts,
      method = rebalance_method,
      freq = move_freq,
      step_days = move_step_days,
      scale_to_time = move_scale_to_time,
      move_mult = move_mult
    )

    if (move_freq == "monthly") {
      ev_times <- which(lubridate::day(dates) == 1)
    } else {
      ev_times <- seq(1, n_days, by = move_step_days)
    }

    out <- ode(
      y = y_init, times = times, func = ode_fun, parms = params,
      method = ode_method, events = list(func = ev_fun, time = ev_times)
    )
  } else {
    out <- ode(
      y = y_init, times = times, func = ode_fun, parms = params,
      method = ode_method
    )
  }

  # Convert to data frame and add date column
  out_df <- as.data.frame(out)

  if (ode_method == "rk4") {
    daily_idx <- match(seq_len(n_days), round(out_df$time, 10))
    if (any(is.na(daily_idx))) {
      stop("RK4 output did not contain all integer daily time points.")
    }
    out_df <- out_df[daily_idx, , drop = FALSE]
    out_df$time <- seq_len(n_days)
  }

  out_df$date <- dates
  out_df
}


#' Run ODE simulation and return only the final state
#'
#' This helper is useful for long burn-ins where we only need the endpoint.
#' To reduce peak memory use, the simulation can be split into chunks and the
#' final state of each chunk passed into the next chunk.
#'
#' @param chunk_size_days Optional chunk size for long runs. If NULL or larger
#'   than the simulation length, the run is done in one block.
#' @return Named state vector from the final time point
run_ode_final_state <- function(MAT_A, MAT_C, dates, y_init, params,
                                MAT_TEMP = NULL,
                                use_movement = TRUE, move_mats = NULL,
                                move_freq = c("monthly", "fixed_step"),
                                move_step_days = 14,
                                move_scale_to_time = TRUE,
                                move_mult = 1,
                                ward_ids, N_star,
                                rebalance_hosts = TRUE,
                                rebalance_method = "scale",
                                use_temp_dependence = FALSE,
                                use_local_coupling = FALSE,
                                W_adj = NULL, adj_deg = NULL,
                                MAT_HABITAT_A = NULL, MAT_HABITAT_C = NULL,
                                chunk_size_days = NULL,
                                ode_method = c("lsoda", "rk4"),
                                rk4_dt = 1) {
  ode_method <- match.arg(ode_method)
  n_days <- length(dates)

  if (is.null(chunk_size_days) || chunk_size_days >= n_days || chunk_size_days < 2) {
    out <- run_ode(
      MAT_A = MAT_A,
      MAT_C = MAT_C,
      MAT_TEMP = MAT_TEMP,
      MAT_HABITAT_A = MAT_HABITAT_A,
      MAT_HABITAT_C = MAT_HABITAT_C,
      dates = dates,
      y_init = y_init,
      params = params,
      use_movement = use_movement,
      move_mats = move_mats,
      move_freq = move_freq,
      move_step_days = move_step_days,
      move_scale_to_time = move_scale_to_time,
      move_mult = move_mult,
      ward_ids = ward_ids,
      N_star = N_star,
      rebalance_hosts = rebalance_hosts,
      rebalance_method = rebalance_method,
      use_temp_dependence = use_temp_dependence,
      use_local_coupling = use_local_coupling,
      W_adj = W_adj,
      adj_deg = adj_deg,
      ode_method = ode_method,
      rk4_dt = rk4_dt
    )
    return(extract_final_state(out))
  }

  current_state <- y_init
  chunk_step <- chunk_size_days - 1
  chunk_starts <- seq(1, n_days - 1, by = chunk_step)

  slice_mat <- function(M, idx) {
    if (is.null(M)) {
      return(NULL)
    }
    M[idx, , drop = FALSE]
  }

  for (chunk_start in chunk_starts) {
    chunk_end <- min(chunk_start + chunk_size_days - 1, n_days)
    idx <- chunk_start:chunk_end

    out_chunk <- run_ode(
      MAT_A = MAT_A[idx, , drop = FALSE],
      MAT_C = MAT_C[idx, , drop = FALSE],
      MAT_TEMP = slice_mat(MAT_TEMP, idx),
      MAT_HABITAT_A = slice_mat(MAT_HABITAT_A, idx),
      MAT_HABITAT_C = slice_mat(MAT_HABITAT_C, idx),
      dates = dates[idx],
      y_init = current_state,
      params = params,
      use_movement = use_movement,
      move_mats = move_mats,
      move_freq = move_freq,
      move_step_days = move_step_days,
      move_scale_to_time = move_scale_to_time,
      move_mult = move_mult,
      ward_ids = ward_ids,
      N_star = N_star,
      rebalance_hosts = rebalance_hosts,
      rebalance_method = rebalance_method,
      use_temp_dependence = use_temp_dependence,
      use_local_coupling = use_local_coupling,
      W_adj = W_adj,
      adj_deg = adj_deg,
      ode_method = ode_method,
      rk4_dt = rk4_dt
    )

    current_state <- extract_final_state(out_chunk)
  }

  current_state
}


# ==============================================================================
# Initial conditions
# ==============================================================================

#' Create initial state vector
#'
#' Sets up endemic-like initial conditions:
#'   - Mostly susceptible hosts with seed infections in random wards
#'   - Low initial vector populations (will grow to equilibrium)
#'
#' @param ward_ids Ward identifiers [n_wards]
#' @param N_star Cattle populations [n_wards]
#' @param params Parameter list
#' @param n_seed_wards Number of wards to seed with infections
#' @param seed_I Number of infectious cattle per seeded ward
#' @return Named vector of initial state [n_wards * 13 compartments]
make_initial_state <- function(ward_ids, N_star, params, n_seed_wards = 10, seed_I = 5) {
  n_wards <- length(ward_ids)

  # Seed infections in random wards
  seed_wards <- sample(ward_ids, min(n_seed_wards, n_wards))
  S0 <- as.numeric(N_star)
  I0 <- rep(0, n_wards)
  I0[ward_ids %in% seed_wards] <- seed_I
  S0 <- pmax(S0 - I0, 0)

  # Initial vector populations (10% of equilibrium, mostly susceptible)
  # At equilibrium: V* = Lambda/mu = lam_max * N * scale / mu
  # Start at 10% of this (will grow during burn-in)
  A_init <- pmax(50, 0.1 * params$lam_A_max * N_star / params$mu_A)
  C_init <- pmax(50, 0.1 * params$lam_C_max * N_star / params$mu_C)

  c(
    # Host compartments
    setNames(S0, nm("Sh", ward_ids)),
    setNames(rep(0, n_wards), nm("Eh", ward_ids)),
    setNames(I0, nm("Ih", ward_ids)),
    setNames(rep(0, n_wards), nm("Rh", ward_ids)),

    # Aedes compartments (99% susceptible, 1% infected to seed transmission)
    setNames(0.99 * A_init, nm("SvA", ward_ids)),
    setNames(rep(0, n_wards), nm("EvA", ward_ids)),
    setNames(0.01 * A_init, nm("IvA", ward_ids)),

    # Culex compartments
    setNames(0.99 * C_init, nm("SvC", ward_ids)),
    setNames(rep(0, n_wards), nm("EvC", ward_ids)),
    setNames(0.01 * C_init, nm("IvC", ward_ids)),

    # Cumulative trackers (start at zero). Hh is split by vector genus so that
    # genus-attributable hazard can be read directly from the ODE output;
    # total host hazard is Hh_A + Hh_C.
    setNames(rep(0, n_wards), nm("Ch",   ward_ids)),
    setNames(rep(0, n_wards), nm("Hh_A", ward_ids)),
    setNames(rep(0, n_wards), nm("Hh_C", ward_ids))
  )
}
