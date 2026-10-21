This repo contains the code to carry out the vector abundance, force of infection and incidence, and dynamic modelling for the manuscript  titled <i>Culex</i>-driven endemicity of Rift Valley fever shapes epidemic risk in northern Tanzania. 

*Jessica Clark1, Tijani A. Sulaimon2, Joshua Carter2, James Mlacha4, Peter Paul3, 
Sarah Cleaveland1, Furaha Mmramba3, Oliva Manangwa3, 
Paul C. D. Johnson1, Jennifer S. Lord2

1.	School of Biodiversity, One Health and Veterinary Medicine, College of Medical, Veterinary and Life Sciences, University of Glasgow, Glasgow, UK
2.	Vector Group, Liverpool School of Tropical Medicine, Pembroke Place, Liverpool, L3 5QA, UK
3.	Vector & Vector-borne Diseases Research Institute, Tanzania Veterinary Laboratory Agency, Tanga, Tanzania

* corresponding author
Jessica.clark@glasgow.ac.uk

The FOI analysis relies on outputs from the vector analysis so the vector script must be run first for both <i>Culex</i> and <i>Aedes</i>. 
The GPS points and covariates have been processed to de-identify the site locations. That means that the outputs from these scripts <b>will not</b> match those within in the manuscript. Please contact the authors if exact GPS points are required, with a short overview of your needs, intentions and a key contact person. These requests will be considered on a case-by-case basis. 

## Dynamic transmission model (`spatial_model/`)

A ward-level, deterministic model of Rift Valley fever virus transmission between
cattle and two mosquito genera (*Aedes* and *Culex*) across 66 wards in Hai, Moshi
and Babati districts, northern Tanzania.

### Model structure

- **Cattle:** SEIRS (susceptible, exposed, infectious, recovered) in each ward, with
  births balancing deaths so herd size stays constant, and waning immunity.
- **Mosquitoes:** SEI for each genus. Recruitment is driven by rainfall: a standing-water
  store for *Aedes* and 14-day rainfall for *Culex*. Biting, mortality and extrinsic
  incubation depend on temperature.
- **Spatial connections:** cattle in neighbouring wards share exposure through local
  coupling (contiguous wards), and cattle trade moves animals between wards monthly.
- **Calibration:** vectors per host (by district and genus) and the *Aedes* water
  half-lives were estimated by approximate Bayesian computation against ward-level
  2023 seroprevalence and 2024 seroconversion. `run_model.R` uses the posterior median
  stored in `output/abc_posterior.rds`.

### Files

| File | Contents |
|---|---|
| `run_model.R` | Runs the model and produces diagnostic output |
| `R/00_config.R` | Parameters, settings and thermal response functions |
| `R/01_data_prep.R` | Loads data and builds the climate, movement and adjacency inputs |
| `R/02_model.R` | Model equations and simulation functions |
| `R/03_helpers.R` | Functions for summarising output and comparing with survey data |
| `climate_data/daily_climate_df.rds` | Daily rainfall and temperature per ward |
| `ward_info/ward_summary_df.rds` | Ward identifiers and cattle population sizes |
| `ward_info/survey_validation.rds` | Ward-level survey targets (seroprevalence and seroconversion) |
| `movement_data/sim_edge_list.rds` | Simulated cattle movements between wards |
| `output/abc_posterior.rds` | Calibrated parameter sets |
| `../SpatialData/contiguous.wards.csv` | Pairs of neighbouring wards |

Survey data are provided as ward-level summaries only. Village-level data are not
included.

### Running the model

Run from inside the `spatial_model/` folder, so the relative file paths resolve:

    cd spatial_model
    Rscript run_model.R

A run takes about 15 minutes on a laptop. It simulates a 40-year burn-in, cycling
2010–2022 climate, to reach an endemic state, then September 2022 to December 2024
with observed climate. It prints summary statistics (prevalence, seroprevalence,
infection probabilities) and saves a validation figure to
`output/fig_validation_median.pdf`.

`run_model.R` uses the posterior median parameter set. Manuscript figures summarise
across all accepted posterior parameter sets, so values from this script will not
match the manuscript exactly.

Required R packages: deSolve, dplyr, tidyr, ggplot2, lubridate, Matrix, tibble, zoo,
patchwork, stringr, scales, igraph, viridis, readr.