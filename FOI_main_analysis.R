#' FOI (2023) and incidence (2023-2024) analyses with main figures, using data_public/.
#' Run genus-level-spatial-analysis-public.R for Culex and Aedes first. Run from the repository root.

rm(list = ls())
library(tidyverse)
library(INLA)
library(fmesher)
library(inlabru)
library(RColorBrewer)
library(lme4)
library(sf)
library(ggspatial)
library(patchwork)
library(colorspace)
library(blockCV)
library(parallel)
require(DescTools)

#' Load local functions
source("functions.R")

dir.create("FOI_analysis/public_outputs", showWarnings = FALSE)

################################################################################
#### Data ####
################################################################################

#' Read cattle serology data including predicted Culex abundance
dat <- read.csv("data_public/cattleSerologyWithCulexPred.csv")
dat <- read.csv("data_public/cattleSerologyWithAedesPred.csv") |>
  dplyr::select(newID, meanAedes) |>
  left_join(dat, by = "newID")

#' Convert geometry column back to sfc
crs.km <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
dat$geometry <- st_as_sfc(dat$geometry, crs = st_crs(crs.km))

#' Convert to sf
dat <- st_as_sf(dat)

#' Read in cattle data with ages
cattle.age <- read.csv("data_public/allCattleDat.csv")
cattle.age$X <- NULL
head(cattle.age)
dat$newAgeInt <- cattle.age$newAgeInt[match(dat$newID, cattle.age$newID)]

#' Both dat and cattle.age have data on serology test results. Check they agree.
dim(dat)
dim(cattle.age)
#' How many cattle in dat are not in cattle.age?
length(setdiff(dat$newID, cattle.age$newID))
#' How many cattle in cattle.age are not in dat?
length(setdiff(cattle.age$newID, dat$newID))
#' Check that test results agree across the two data sets
table(cattle.age$status2023[match(dat$newID, cattle.age$newID)],
      dat$exposed)
#' What about the animals that were followed up in 2024? How many
#' and what were the results?
table(cattle.age$status2024[match(dat$newID, cattle.age$newID)],
      exclude = NULL)

dat$exposed24 <- as.integer(cattle.age$status2024 == "Positive")[match(dat$newID, cattle.age$newID)]
dat$exposed24[dat$exposed == 1] <- NA
dat$exposed.yrs <- 8/12 # followed up after 8 months

colSums(is.na(dat))

dat$culexAbundance.sds <- scale(dat$meanCulex)
dat$aedesAbundance.sds <- scale(dat$meanAedes)

#' Make factors
dat$District <- factor(dat$District)
dat$Village <- factor(dat$Village)
dat$hhID <- factor(dat$hhID)
dat$Class <- factor(dat$Class, c("Rural", "Periurban"))
dat$ClassPeriurban <- as.integer(dat$Class == "Periurban")

# sort out villages
cattle.age$Village[which(cattle.age$Village=="Gidewali")] <- "Gidawari"

dat <- st_as_sf(dat)

################################################################################
#### Modelling: 2023 force of infection ####
################################################################################

#' Store the model inputs and outputs in a list
mods <-
  list(priors = list(),
       inputs = list(),
       outputs = list())

#' Add the priors to the list. These will be the same for all models.
#' Prior on IID random effect variances.
mods$priors$iid.prec.prior <- list(prec = list(prior = "loggamma", param = c(0.5, 0.01)))
#' Prior on the intercept effect. This is the same as the default prior on other
#' fixed effects.
mods$priors$prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.001)


#### Spatial model bits ####

#' Set up the mesh
set.seed(123) # so that we always get the same mesh
mods$inputs$spatial$hull <- fm_extensions(dat, convex = c(25, 50))
mods$inputs$spatial$mesh <- fm_mesh_2d(dat, max.edge = c(7, 40), cutoff = 0.005,
                                       boundary = mods$inputs$spatial$hull)

#' Re-set random seed so that we notice if the model results are stochastic
rand.seed <- round((as.numeric(Sys.time()) - floor(as.numeric(Sys.time())))*10000000)
set.seed(rand.seed)

#' inla.spde2.matern() makes a Matern SPDE model object
#' The SPDE allows us to approximate in a computationally efficient way
#' a continuous Matern correlation structure with a mesh.
mods$inputs$spatial$spde <-
  inla.spde2.pcmatern(mesh = mods$inputs$spatial$mesh,
                      prior.range = c(1, 0.025),
                      prior.sigma = c(1, 0.1),
                      constr = TRUE)

#' Specify the spatial model components
mods$inputs$spatial$components <-
  ~ 0 + Intercept(1) +
  newAgeInt(newAgeInt, model = "offset") +
  culexAbundance.sds(culexAbundance.sds, model = "linear")  +
  aedesAbundance.sds(aedesAbundance.sds, model = "linear")  +
  ndvi.sds(ndvi.sds, model = "linear")  +
  ndvi.sds2(ndvi.sds2, model = "linear")  +
  NRainDays.log10.sds(NRainDays.log10.sds, model = "linear")  +
  proportionIrrigatedCrop.sds(proportionIrrigatedCrop.sds, model = "linear")  +
  elevation500m.sds(elevation500m.sds, model = "linear")  +
  ClassPeriurban(ClassPeriurban, model = "linear")  +
  popDens500m.log10.sds(popDens500m.log10.sds, model = "linear")  +
  buildingDensity500m.log10.sds(buildingDensity500m.log10.sds, model = "linear")  +
  Village(Village, model = "iid", hyper = mods$priors$iid.prec.prior) +
  hhID(hhID, model = "iid", hyper = mods$priors$iid.prec.prior) +
  field(geometry, model = mods$inputs$spatial$spde)


#### Null spatial fixed effects model ####
mods$inputs$spatial.null <- mods$inputs$spatial
mods$inputs$spatial.null$obs <-
  bru_obs(
    formula = exposed ~ Intercept + Village + hhID + field + offset(log(newAgeInt)),
    family = "binomial",
    control.family = list(link = "cloglog"),
    data = dat)

#' fit the model
mods$outputs$spatial.null <-
  bru(
    components = mods$inputs$spatial.null$components,
    obs = mods$inputs$spatial.null$obs,
    options = list(
      control.fixed = mods$priors$prior.fixed,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                             config = TRUE,
                             return.marginals.predictor = TRUE)))

summary(mods$outputs$spatial.null)


# plot the spatial range what range does the autocorrelation fade?
autocorrealtionPlot <-
  ggplot(data = spde.posterior(mods$outputs$spatial.null, "field", what = "matern.correlation")) +
  geom_ribbon(aes(x = x, ymin = `q0.025`, ymax = `q0.975`), fill = "grey", alpha = 0.2) +
  geom_line(aes(x = x, y = q0.5)) +
  xlab("Distance (km)") +
  ylab("Matern correlation") +
  scale_x_continuous(breaks = seq(0,60, by = 5)) +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(2, 2, 2, 2)) +
  geom_hline(yintercept = 0.14, lty = 2) +
  geom_vline(xintercept = mods$outputs$spatial.null$summary.hyperpar["Range for field", "0.5quant"], lty = 2)
autocorrealtionPlot

#### Spatial GLMM (w/ vectors)####
mods$inputs$spatialVectors <- mods$inputs$spatial
#' Specify the observation model
mods$inputs$spatialVectors$obs <-
  bru_obs(
    formula = exposed ~ Intercept +
      culexAbundance.sds +
      aedesAbundance.sds +
      Village +
      hhID +
      field +
      offset(log(newAgeInt)),
    family = "binomial",
    control.family = list(link = "cloglog"),
    data = dat)

#' Fit the model
mods$outputs$spatialVectors <-
  bru(
    components = mods$inputs$spatialVectors$components,
    obs = mods$inputs$spatialVectors$obs,
    options = list(
      control.fixed = mods$priors$prior.fixed,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                             config = TRUE,
                             return.marginals.predictor = TRUE)))
summary(mods$outputs$spatialVectors)

# plot to assess
fit.summary.fixed <- mods$outputs$spatialVectors$summary.fixed
fit.summary.fixed$ID <- rownames(fit.summary.fixed)
fit.summary.fixed$ID <- factor(fit.summary.fixed$ID, fit.summary.fixed$ID)

vectorFOIspatial <-
  fit.summary.fixed %>%
   mutate(ID = fct_recode(ID
                          ,Culex = "culexAbundance.sds"
                          ,Aedes = "aedesAbundance.sds"
   )) %>%
  ggplot(aes(y = ID, x = exp(mean))) +
  geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
  geom_vline(xintercept = 1, linetype = 2, colour = "lightgrey") +
   scale_y_discrete(labels = c(
    "Intercept"
    ,expression(italic("Culex"))
    ,expression(italic("Aedes"))
   )) +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(2, 2, 2, 2)) +
  labs(y = "", x = "Estimate ± 95% CI"
       )
vectorFOIspatial

##### All environmental covars associated with abundance ####

mods$inputs$allSpatialCovars <- mods$inputs$spatial
mods$inputs$allSpatialCovars$obs <-
  bru_obs(
    formula =
      exposed ~
      Intercept +
      ndvi.sds +
      ndvi.sds2 +
      NRainDays.log10.sds +
      elevation500m.sds +
      buildingDensity500m.log10.sds +
      popDens500m.log10.sds +
      proportionIrrigatedCrop.sds +
      Village +
      hhID +
      field +
      offset(log(newAgeInt)),
    family = "binomial",
    control.family = list(link = "cloglog"),
    data = dat)


# fit the model
mods$outputs$allSpatialCovars <-
  bru(
    components = mods$inputs$allSpatialCovars$components,
    obs = mods$inputs$allSpatialCovars$obs,
    options = list(
      control.fixed = mods$priors$prior.fixed,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                             config = TRUE,
                             return.marginals.predictor = TRUE)))

summary(mods$outputs$allSpatialCovars)

# visualise this
fit.summary.fixed <- mods$outputs$allSpatialCovars$summary.fixed
fit.summary.fixed$ID <- rownames(fit.summary.fixed)
fit.summary.fixed$ID <- factor(fit.summary.fixed$ID, fit.summary.fixed$ID)
fit.summary.fixed <-
  fit.summary.fixed %>%
  mutate(ID = fct_recode(ID,
                         NDVI = "ndvi.sds",
                         `NDVI²` = "ndvi.sds2",
                         `proportion \n irrigated` = "proportionIrrigatedCrop.sds",
                         elevation = "elevation500m.sds",
                         `population density` = "popDens500m.log10.sds",
                         `building density` = "buildingDensity500m.log10.sds",
                         `n rainy days` = "NRainDays.log10.sds"))

# Plot
all.spatial.covars.plot <-
  fit.summary.fixed %>%
  ggplot(aes(y = ID, x = exp(mean))) +
  geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
  geom_vline(xintercept = 1, linetype = 2) +
  theme_minimal(base_size = 14
                ) +
  labs(y = "", x = "Estimate ± 95% CI"
       )

all.spatial.covars.plot
sapply(mods$outputs, function(f) f$waic$waic)

### Blocked k-fold cross-validation ####

#' We are using blocked cross-validation due to the non-independence problem
#' with randomly selecting held-out folds:
#' Roberts et al. (2017). Cross-validation strategies for data with temporal,
#' spatial, hierarchical, or phylogenetic structure. Ecography.
#' Valavi et al. (2019). blockCV: An R package for generating spatially or
#' environmentally separated folds for cross-validation of species distribution
#' models. Methods in Ecology and Evolution.

#' The aim is to choose clusters that are spatially and temporally independent.
#' We do this by first estimating the spatial scale of correlation
mods$outputs$spatial.null$summary.hyperpar["Range for field", "mean"]
#' ~10 km
#' A robust way to do this is by k-means clustering.
#' We want the diameter of the held-out clusters of livestock sampling to be at least double that scale.
#' More folds (higher k) will give smaller clusters. 5 and 10 are common choices for k.
#' I tried both, and only k=5 gave a minimum cluster diameter that was > 4 * 10 km. The code
#' for calculating cluster diameter (2 * distance from fold point centroid to nearest
#' outside point) is in the CV loop below.

#' Assign clusters (the fold maps are in the supplementary script)
cluster.k <- 5
set.seed(456)
cv <- cv_cluster(x = dat, k = cluster.k)

#' Store cluster ID in dat.all for now (don't want to do anything silly to the main df)
dat.all <- dat
dat.all$fold <- NA
dat.all$fold[1:nrow(cv$biomod_table)] <-
  ((1 - cv$biomod_table) %*% 1:cluster.k)[, 1]
table(dat.all$fold)

#' Set inla options for CV

inla.opt <-
  list(
    num.threads = "1:1",
    inla.mode = "experimental",
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                           config = TRUE,
                           return.marginals.predictor = TRUE,
                           openmp.strategy = NULL),
    control.family = list(link = "cloglog"))

inla.opt.cv <- inla.opt
inla.opt.cv$control.compute$cpo <- inla.opt.cv$control.compute$dic <-
  inla.opt.cv$control.compute$waic <- FALSE
inla.opt.cv$control.compute$return.marginals.predictor <- FALSE
inla.opt.cv$num.threads <- detectCores()
inla.opt.cv$inla.mode <- "compact"
inla.opt.cv$control.inla = list(strategy = "simplified.laplace", int.strategy = "eb", h = 0.02)


modelSpecs <-
  list(spatial.null =  "offset(log(newAgeInt))")
modelSpecs$spatialVectors <-  "culexAbundance.sds + aedesAbundance.sds + offset(log(newAgeInt))"
modelSpecs$allSpatialCovars <-
   "ndvi.sds + ndvi.sds2 + NRainDays.log10.sds + elevation500m.sds + buildingDensity500m.log10.sds + popDens500m.log10.sds + proportionIrrigatedCrop.sds  + offset(log(newAgeInt))"


#' Loop over each fold, fitting model to training folds and predicting
#' FOI on held-out fold

#' Storage for the CV diagnostic checks below: per-model out-of-fold
#' predictions and per-fold Brier scores, so we can check (1) whether models
#' actually make different predictions and (2) whether the Brier score
#' ranking is consistent across folds rather than driven by one or two.
cv.diagnostics.2023 <- list()

cv.results.2023 <-
  sapply(names(mods$outputs)[1:3], function(fit.name) {

    dat.all$pred.fold <- NA
    dat.all$fold.diameter <- NA

    print(fit.name)
    fit <- mods$outputs[[fit.name]]

    #' Original components
    components.cv <- mods$inputs[[fit.name]]$components

    #' Loop over all k clusters
    for(i in 1:cluster.k) {
      print(i)
      print(table(!is.na(dat.all$pred.fold)))

      # Calculate distance from centroid of fold to nearest point outside fold
       # get approximate diameter for the distance from the fold centroid to the nearest point outside the fold (scale of fold)
      fold.centroid <- colMeans(st_coordinates(dat.all[dat.all$fold %in% i, ]))
      dat.all$fold.diameter[dat.all$fold %in% i] <-
        2 * min(st_distance(st_sfc(st_point(fold.centroid), crs = crs.km),
                            dat.all[dat.all$fold %in% (1:ncol(cv$biomod_table))[-i], ]))
      print(unique(dat.all$fold.diameter[dat.all$fold %in% i]))

      # Refit model to 4 of the 5 folds, leaving 1 fold out (fold i)
      inla.opt.cv$verbose <- inla.opt$verbose <- FALSE
      dat.all.cv <- dat.all
      dat.all.cv$exposed[dat.all.cv$fold %in% i] <- NA

      # Specify the observation model
      obs.cv <-
        bru_obs(
          formula = formula(paste0("exposed ~ Intercept + Village + hhID + field + ", modelSpecs[which(names(modelSpecs)==fit.name)])),
          family = "binomial",
          data = dat.all.cv)

      # Fit the model
      cv.fit.start <- Sys.time()
      fit.cv <-
        bru(
          components = components.cv,
          obs = obs.cv,
          options = inla.opt)
      cv.fit.finish <- Sys.time()
      print(cv.fit.finish - cv.fit.start)

      # Predict left-out fold from model fitted to other 4 folds.
      dat.all$pred.fold[which(!cv$biomod_table[, i])] <-
        predict(fit.cv,
                formula = paste0("Intercept + ", modelSpecs[[fit.name]]),
                newdata = dat.all.cv)$mean[which(!cv$biomod_table[, i])]
        print(table(!is.na(dat.all$pred.fold)))
    }

    #' Back-transform from the cloglog-scale linear predictor to a probability
    dat.all$probPred <- 1 - exp(-exp(dat.all$pred.fold))
    plot(dat.all$probPred, dat.all$exposed)

    brierScore <- BrierScore(dat.all$exposed, pred = dat.all$probPred)

    #' Brier score per fold, so we can check afterwards whether the overall
    #' ranking is consistent across folds or driven by just one or two
    foldBrier <-
      sapply(1:cluster.k, function(k) {
        rows <- which(dat.all$fold == k)
        BrierScore(dat.all$exposed[rows], pred = dat.all$probPred[rows])
      })
    cv.diagnostics.2023[[fit.name]] <<- list(probPred = dat.all$probPred, foldBrier = foldBrier)

    brierScore
  })

##### CV diagnostic checks #####

#' 1. Do the models' out-of-fold predicted probabilities actually differ?
#' If they're near-identical, a lower Brier score for the simpler model isn't
#' meaningful - it just means the extra covariates aren't changing anything.
predMatrix <- sapply(cv.diagnostics.2023, function(x) x$probPred)
cat("Correlation between models' out-of-fold predicted probabilities:\n")
print(cor(predMatrix, use = "pairwise.complete.obs"))

write.csv(print(cor(predMatrix, use = "pairwise.complete.obs")), "FOI_analysis/public_outputs/predCorrMatrix.csv")
#' 2. Is the Brier score ranking consistent across all folds, or is it
#' driven by one or two folds?
foldBrierTab <- sapply(cv.diagnostics.2023, function(x) x$foldBrier)
rownames(foldBrierTab) <- paste0("fold", 1:cluster.k)
cat("\nBrier score by fold and model:\n")
print(round(foldBrierTab, 3))
write.csv(foldBrierTab, "FOI_analysis/public_outputs/foldBScores.csv")


cat("\nHow many folds does each model 'win' (lowest Brier)?\n")
foldWinners <- factor(apply(foldBrierTab, 1, function(x) names(x)[which.min(x)]),
                      levels = colnames(foldBrierTab))
print(table(foldWinners))

write.csv(foldWinners, "FOI_analysis/public_outputs/foldWins.csv")

#' ## Tabulate results for all models

#' Get WAIC for all models

WAIC.tab <- cbind(Fixed_effects = unlist(modelSpecs), IC.inla(mods$outputs[1:3], crit = "waic"))

covariates <-
  c(NDVI = "ndvi",
    NDVI2 = "ndvi.sds2",
    NRainDays = "NRainDays.log10",
    PropIrrigated = "proportionIrrigatedCrop",
    Elevation = "elevation500m",
    BuildDens = "buildingDensity500m.log10",
    PopDens = "popDens500m.log10",
    Culex = "culexAbundance.sds",
    Aedes = "abundanceAedes.sds" )

for(i in 1:nrow(WAIC.tab)) {
  for(j in 1:length(covariates)) {
    WAIC.tab[i, "Fixed_effects"] <-
      gsub(covariates[j], names(covariates)[j], WAIC.tab[i, "Fixed_effects"])
  }
  WAIC.tab[i, "Fixed_effects"] <- gsub("\\.sds", "", WAIC.tab[i, "Fixed_effects"])
}
names(WAIC.tab) <- gsub("delta", "Δ", names(WAIC.tab))

#' Add 5-fold CV results to WAIC table

#The lower the Brier score is for a set of predictions, the better the predictions are calibrated.
WAIC.tab <- cbind(WAIC.tab, as.matrix(round(cv.results.2023, 4))[rownames(WAIC.tab), ])
colnames(WAIC.tab) <- c("Fixed_effects"
                        ,"WAIC"
                        ,"ΔWAIC"
                        ,"Weights"
                        ,"BrierScore")

write.csv(WAIC.tab,  "FOI_analysis/public_outputs/modelFitTab.csv")

################################################################################
#### Extract the posteriors to get Village-level FOI ####
################################################################################

X <- as.data.frame(model.matrix(~ newAgeInt +
                                  ndvi.sds +
                                  ndvi.sds2 +
                                  proportionIrrigatedCrop.sds +
                                  elevation500m.sds +
                                  buildingDensity500m.log10.sds +
                                  popDens500m.log10.sds +
                                  NRainDays.log10.sds,
                                data = dat))

#' Which estimate to use? Mean or median? Try each and
#' assess sensitivity to choice.

lower.est <- "0.025quant"
upper.est <- "0.975quant"

point.est <- c("mean", "0.5quant")[2]

#' Evaluate the spatial field (object field.est) from our null model only at our sampled locations
spatialRE <-
  fm_evaluate(fm_evaluator(mods$inputs$spatial.null$mesh, loc = dat),
              mods$outputs$spatial.null$summary.random$field[, point.est])

#' Sum the fixed effects predictions, the Village random effect, and the spatial
#' random effect, giving mean Village-level log FOI predictions
fixedPred <-
  (as.matrix(X[, 1]) %*% mods$outputs$spatial.null$summary.fixed[, point.est])

villageRE <-
  mods$outputs$spatial.null$summary.random$Village[
    match(dat$Village, mods$outputs$spatial.null$summary.random$Village$ID), point.est]

dat$logvillageFOI <-
  fixedPred +
  villageRE +
  spatialRE

#' Sum the lower bound estimates for the fixed effects predictions,  Village random effect, and  spatial
#' random effect, giving 2.5% CI Village-level log FOI predictions
spatialRELow <-
  fm_evaluate(fm_evaluator(mods$inputs$spatial.null$mesh, loc = dat),
              mods$outputs$spatial.null$summary.random$field[, lower.est])
fixedPredLow <-
  (as.matrix(X[, 1]) %*% mods$outputs$spatial.null$summary.fixed[, lower.est])
villageRELow <-
  mods$outputs$spatial.null$summary.random$Village[
    match(dat$Village, mods$outputs$spatial.null$summary.random$Village$ID), lower.est]

dat$logvillageFOILow <-
  fixedPredLow +
  villageRELow +
  spatialRELow

#' Sum the upper bound estimates for the fixed effects predictions,  Village random effect, and  spatial
#' random effect, giving 97.5% CI Village-level log FOI predictions
spatialREHigh <-
  fm_evaluate(fm_evaluator(mods$inputs$spatial.null$mesh, loc = dat),
              mods$outputs$spatial.null$summary.random$field[, upper.est])
fixedPredHigh <-
  (as.matrix(X[, 1]) %*% mods$outputs$spatial.null$summary.fixed[, upper.est])

villageREHigh <-
  mods$outputs$spatial.null$summary.random$Village[
    match(dat$Village, mods$outputs$spatial.null$summary.random$Village$ID), upper.est]

dat$logvillageFOIHigh <-
  fixedPredHigh +
  villageREHigh +
  spatialREHigh

#' Add the herd random effect, giving herd-level log FOI predictions
herdRE <-
  mods$outputs$spatial.null$summary.random$hhID$`0.5quant`[
    match(as.character(dat$hhID), mods$outputs$spatial.null$summary.random$hhID$ID)]
dat$logherdFOI <- dat$logvillageFOI + herdRE


#' Exponentiate log FOI to FOI
dat$villageFOI <- exp(dat$logvillageFOI)
dat$herdFOI <- exp(dat$logherdFOI)

dat$villageFOIlower <- exp(dat$logvillageFOILow)
dat$villageFOIhigher <- exp(dat$logvillageFOIHigh)

#' This is the RAW foi. What this is basically saying, is an animal would experienced
#' FOI number of of infectious exposure events per year if it kept being exposed indefinitely
#' but this isn't possible because seroconverting is a one time event (in this framewokr)
#' so basically how many exposure events happen on average...

data.frame(
  quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
  villageFOI = quantile(dat$villageFOI, probs = c(0, 0.025, 0.5, 0.975, 1)),
  villageFOIlower = quantile(dat$villageFOIlower, probs = c(0, 0.025, 0.5, 0.975, 1)),
  villageFOIhigher = quantile(dat$villageFOIhigher, probs = c(0, 0.025, 0.5, 0.975, 1)),
  row.names = NULL)

#' But what we want to know actually is what is the chance at least one
#' event happens before the year is out?
#' so this is where this comes in....

#' Annual % of susceptible animals expected to seroconvert (bounded 0-100%,
#' unlike the raw FOI rate above): 1 - exp(-FOI), expressed as a percentage.
#' although it is actually the probability of seroconverting per year.
#' 1 - exp(-x) is monotonically increasing, so the lower/upper FOI bounds map
#' to the lower/upper bounds of the percentage without swapping.
dat$annualPercentSeroconversion <- (1-exp(-dat$villageFOI))*100
dat$annualPercentSeroconversionLower <- (1 - exp(-dat$villageFOIlower)) * 100
dat$annualPercentSeroconversionHigher <- (1 - exp(-dat$villageFOIhigher)) * 100

data.frame(
  quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
  annualvillageSC = quantile(dat$annualPercentSeroconversion, probs = c(0, 0.025, 0.5, 0.975, 1)),
  annualvillageScLow = quantile(dat$annualPercentSeroconversionLower, probs = c(0, 0.025, 0.5, 0.975, 1)),
  annualvillageScHigh = quantile(dat$annualPercentSeroconversionHigher, probs = c(0, 0.025, 0.5, 0.975, 1)),
  row.names = NULL)

#' Same summary, broken down by District
dat %>%
  st_drop_geometry() %>%
  group_by(District) %>%
  summarise(
    n = n(),
    `0%` = quantile(annualPercentSeroconversion, 0),
    `2.5%` = quantile(annualPercentSeroconversion, 0.025),
    `50%` = quantile(annualPercentSeroconversion, 0.5),
    `97.5%` = quantile(annualPercentSeroconversion, 0.975),
    `100%` = quantile(annualPercentSeroconversion, 1)
  )

#' at low FOIs FOI ~= to the annual percent, but they diverge at higher FOIs.
#' so an expected count (i.e. incident rate)

#' Express FOI as a standard Poisson incidence rate (cases per 100 animal-years
#' at risk). FOI is already a per-animal-year rate because age (in years) was
#' used as the model offset, so this is a re-expression, not a new estimate:
#' N infections in T animal-years of exposure ~ Poisson(FOI * T)
incidenceRate.tab <-
  data.frame(
    quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
    incidenceRate100AY = quantile(dat$villageFOI, probs = c(0, 0.025, 0.5, 0.975, 1)) * 100,
    row.names = NULL)

write.csv(incidenceRate.tab, "FOI_analysis/public_outputs/incidenceRatePer100AnimalYears.csv", row.names = FALSE)

#' It's clear that the herd random effect makes very little
#' difference to the Village FOI. Not surprising give the
#' high precision (low variance) estimates. This is also reflected
#' in the variances of the model components:
var(fixedPred)

#' ...giving a high marginal R-squared on the log scale (ignoring binomial variation)
var(fixedPred) / sum(var(fixedPred), var(villageRE), var(herdRE), var(spatialRE))

#' Calculate each animal's predicted probability of seropositivity
#' using probPos <- 1-exp(-exp(log(age) + log(FOI)))
dat$probPos <- 1-exp(-exp(log(dat$newAgeInt) + log(dat$herdFOI)))

###### Village-level predictions ######

VillagePred <-
  data.frame(Village = sort(unique(dat$Village)))
rownames(VillagePred) <- VillagePred$Village
VillagePred[, c("X", "Y")] <-
  do.call("rbind",
          tapply(dat$geometry, dat$Village,
                 function(x) apply(st_coordinates(x), 2, mean))[VillagePred$Village])
VillagePred <- st_as_sf(VillagePred, coords = c("X", "Y"), crs = st_crs(dat), remove = FALSE)
VillagePred$District <- factor(tapply(as.character(dat$District), dat$Village, unique)[as.character(VillagePred$Village)])
VillagePred$Village <- factor(VillagePred$Village)
table(VillagePred$District)


VillagePred$NTested <-
  as.vector(table(dat$Village)[VillagePred$Village])
VillagePred$obsPos <-
  tapply(dat$exposed, dat$Village, sum)[VillagePred$Village]
VillagePred$predPos <-
  tapply(dat$probPos, dat$Village, sum)[VillagePred$Village]
VillagePred$FOI <-
  tapply(dat$villageFOI, dat$Village, mean)[VillagePred$Village]

#' Annual % of susceptible animals expected to seroconvert (bounded 0-100%,
#' unlike the raw FOI rate above): 1 - exp(-FOI), expressed as a percentage
VillagePred$annualPercentSeroconversion <-
  (1 - exp(-VillagePred$FOI)) * 100

VillagePred$meanAge <-
  tapply(dat$newAgeInt, dat$Village, mean)[VillagePred$Village]

VillagePred$culexAbundance <-
  tapply(dat$culexAbundance, dat$Village, mean)[VillagePred$Village]

VillagePred$aedesAbundance <-
  tapply(dat$aedesAbundance, dat$Village, mean)[VillagePred$Village]

VillagePred <- VillagePred %>%
  left_join(
    dat %>%
      st_drop_geometry() %>%
      dplyr::select(Village, annualPercentSeroconversionLower,
                    annualPercentSeroconversionHigher) %>%
      group_by(Village) %>%
      summarise(meanAnnLow = mean(annualPercentSeroconversionLower),
                meanAnnHigh = mean(annualPercentSeroconversionHigher)),
    by = "Village")

###### Make FOI plot ######

(FOIPlot <-
  ggplot() +
  annotation_map_tile(type = "cartolight", zoomin = 0) +  # You can change the tile type
  geom_sf(data = VillagePred,
          aes(fill = annualPercentSeroconversion, shape = District),
          size = 5, colour = "darkgrey") +
  annotation_scale(location = "tl", width_hint = 0.5) +
  coord_sf(crs = crs.km, expand = T) +
  scale_fill_viridis_c(option = "A",
                        direction = -1,
                        breaks = scales::breaks_width(2)) +
  scale_shape_manual(values = c(Babati = 21, Hai = 22, Moshi = 24)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "bottom",
    legend.title = element_text(size = 12)
  )+
  labs(fill = "Annual % of \nsusceptible \nseroconversions") )


data.frame(
  quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
  medianAnnVill = quantile(VillagePred$annualPercentSeroconversion, probs = c(0, 0.025, 0.5, 0.975, 1)),
  row.names = NULL)

# write csv to send to Tijani
write.csv(VillagePred %>%
            select(-geometry), "FOI_analysis/public_outputs/VillagePredictions.csv", row.names = F)

##### Main figure: FOI panel #####
Culex_pred_count_mean <- readRDS("figures_public/Culex_pred_count_mean.rds")
Aedes_pred_count_mean <- readRDS("figures_public/Aedes_pred_count_mean.rds")

Culex_pred_count_mean[[1]]@labels$title <- ""
Aedes_pred_count_mean[[1]]@labels$title <- ""

#' Culex and Aedes maps already share the same scale_colour_gradientn()
#' limits/breaks/palette, so one legend correctly represents both. Pull it
#' from Culex's main panel ([[1]]; [[2]] is the inset) and drop it from both.
sharedVectorLegend <- cowplot::get_legend(
  Culex_pred_count_mean[[1]] +
    guides(colour = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(0.4, "cm"))) +
    theme(legend.position = "bottom",
          plot.margin = margin(2, 2, 2, 2),
          legend.title = element_text(size = 12))
)

Culex_noleg <- Culex_pred_count_mean
Culex_noleg[[1]] <- Culex_noleg[[1]] + theme(legend.position = "none")
Aedes_noleg <- Aedes_pred_count_mean
Aedes_noleg[[1]] <- Aedes_noleg[[1]] + theme(legend.position = "none")

topRow <- cowplot::plot_grid(
  cowplot::plot_grid(Culex_noleg, Aedes_noleg, ncol = 2, labels = c("A", "B")),
  sharedVectorLegend,
  nrow = 2, rel_heights = c(1, 0.15)
)

bottomGrid <- cowplot::plot_grid(
  FOIPlot, autocorrealtionPlot, vectorFOIspatial, all.spatial.covars.plot,
  nrow = 2, ncol = 2, labels = c("C", "D", "E", "F")
)

cowplot::plot_grid(topRow, bottomGrid, nrow = 2, rel_heights = c(1, 2))
ggsave("FOI_analysis/public_outputs/FOIpanel.pdf")

################################################################################
#### Incidence analysis (2023-2024 seroconversion) ####
################################################################################

# load incidence data

cattleDatInc <- read.csv("data_public/cattleIncidenceDataWithAedesPred.csv")
cattleDatInc <- read.csv("data_public/cattleIncidenceDataWithCulexPred.csv") |>
  dplyr::select(newID, meanCulex) %>%
  left_join(cattleDatInc, by = "newID") |>
  mutate(culexAbundance.sds = scale(meanCulex),
    aedesAbundance.sds = scale(meanAedes))

cattleDatInc <- st_as_sf(cattleDatInc, coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(crs.km)

#' District isn't in the incidence CSVs, so bring it in from cattle.age
#' (one District per Village, so this is a safe 1:1 join -- verified no
#' Village maps to >1 District and no animal is left unmatched)
cattleDatInc <- cattleDatInc %>%
  left_join(cattle.age %>% distinct(Village, District), by = "Village")

#'Fit serocatalytic model to incidence data for seronegative animals re-sampled in 2024
#'Does vector abundance explain the 2024 seroconversion in light of El Nino

mods$inputs$spatial24 <- mods$inputs$spatial
mods$inputs$spatial24$obs <-
  bru_obs(
    formula =
      exposed ~
      Intercept +
      culexAbundance.sds +
      aedesAbundance.sds +
      Village +
      hhID +
      field +
      offset(log(prop_year)),
    family = "binomial",
    control.family = list(link = "cloglog"),
    data = cattleDatInc)


#' Fit model
mods$outputs$spatial24 <-
  bru(
    components = mods$inputs$spatial24$components,
    obs = mods$inputs$spatial24$obs,
    options = list(
      control.fixed = mods$priors$prior.fixed,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                             config = TRUE,
                             return.marginals.predictor = TRUE)))
summary(mods$outputs$spatial24)
exp(mods$outputs$spatial24$summary.fixed)

# visualise this
fit.summary.fixed <- mods$outputs$spatial24$summary.fixed
fit.summary.fixed$ID <- rownames(fit.summary.fixed)
fit.summary.fixed$ID <- factor(fit.summary.fixed$ID, fit.summary.fixed$ID)

# Plot
vectorIncidence <-
  fit.summary.fixed %>%
  mutate(ID = fct_recode(ID
                         ,Culex = "culexAbundance.sds"
                         ,Aedes = "aedesAbundance.sds"
  )) %>%
  ggplot(aes(y = ID, x = exp(mean))) +
  geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
  geom_vline(xintercept = 1, linetype = 2) +
   scale_y_discrete(labels = c(
     "Intercept"
     ,expression(italic("Culex"))
     ,expression(italic("Aedes"))
   )) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")+
  labs(y = "", x = "Estimate ± 95% CI")

###### Null spatial fixed effects model #####
mods$inputs$spatial24.null <- mods$inputs$spatial
mods$inputs$spatial24.null$obs <-
  bru_obs(
    formula =
      exposed ~
      Intercept +
      Village + hhID + field + offset(log(prop_year)),
    family = "binomial",
    control.family = list(link = "cloglog"),
    data = cattleDatInc)


#' Fit spatial model
mods$outputs$spatial24.null <-
  bru(
    components = mods$inputs$spatial24.null$components,
    obs = mods$inputs$spatial24.null$obs,
    options = list(
      control.fixed = mods$priors$prior.fixed,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                             config = TRUE,
                             return.marginals.predictor = TRUE)))

summary(mods$outputs$spatial24.null)

exp(mods$outputs$spatial24.null$summary.fixed)

sapply(mods$outputs, function(f) f$waic$waic)


#####  Apply CV to incidence modelling  #####

#' Assign clusters
#' cv_cluster() uses k-means internally, so fix the seed - otherwise the fold
#' assignment (and therefore the CV comparison between models) is different
#' every time this is re-run
cluster.k <- 5
set.seed(456)
cv <- cv_cluster(x = dat, k = cluster.k)

#' Store cluster ID in dat.all for now (don't want to do anything silly to the main df)
dat.all.inc <- dat
dat.all.inc$fold <- NA
dat.all.inc$fold[1:nrow(cv$biomod_table)] <-
  ((1 - cv$biomod_table) %*% 1:cluster.k)[, 1]
table(dat.all.inc$fold)

inla.opt <-
  list(
    num.threads = "1:1",
    inla.mode = "experimental",
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE,
                           config = TRUE,
                           return.marginals.predictor = TRUE,
                           openmp.strategy = NULL),
    control.family = list(link = "cloglog"))

inla.opt.cv <- inla.opt
inla.opt.cv$control.compute$cpo <- inla.opt.cv$control.compute$dic <-
  inla.opt.cv$control.compute$waic <- FALSE
inla.opt.cv$control.compute$return.marginals.predictor <- FALSE
inla.opt.cv$num.threads <- detectCores()
inla.opt.cv$inla.mode <- "compact"
inla.opt.cv$control.inla = list(strategy = "simplified.laplace", int.strategy = "eb", h = 0.02)


modelSpecsIncidence <-
  list(spatial24 = "culexAbundance.sds + aedesAbundance.sds + offset(log(exposed.yrs))")
modelSpecsIncidence$spatial24.null <-
  "offset(log(exposed.yrs))"


#' Loop over each fold, fitting model to training folds and predicting
#' FOI on held-out fold

#' Storage for the CV diagnostic checks below: per-model out-of-fold
#' predictions and per-fold Brier scores, so we can check (1) whether models
#' actually make different predictions and (2) whether the Brier score
#' ranking is consistent across folds rather than driven by one or two.
cv.diagnostics <- list()

cv.results <-
  sapply(names(mods$outputs)[4:length(names(mods$outputs))], function(fit.name) {

    dat.all.inc$pred.fold <- NA
    dat.all.inc$fold.diameter <- NA

    print(fit.name)
    fit <- mods$outputs[[fit.name]]

    #' Original components
    components.cv <- mods$inputs[[fit.name]]$components

    #' Loop over all k clusters
    for(i in 1:cluster.k) {
      print(i)
      print(table(!is.na(dat.all.inc$pred.fold)))

      # Calculate distance from centroid of fold to nearest point outside fold
      # get approximate diameter for the distance from the fold centroid to the nearest point outside the fold (scale of fold)
      fold.centroid <- colMeans(st_coordinates(dat.all.inc[dat.all.inc$fold %in% i, ]))
      dat.all.inc$fold.diameter[dat.all.inc$fold %in% i] <-
        2 * min(st_distance(st_sfc(st_point(fold.centroid), crs = crs.km),
                            dat.all.inc[dat.all.inc$fold %in% (1:ncol(cv$biomod_table))[-i], ]))
      print(unique(dat.all.inc$fold.diameter[dat.all.inc$fold %in% i]))

      # Refit model to 4 of the 5 folds, leaving 1 fold out (fold i)
      inla.opt.cv$verbose <- inla.opt$verbose <- FALSE
      dat.all.cv <- dat.all.inc
      dat.all.cv$exposed24[dat.all.cv$fold %in% i] <- NA

      # Specify the observation model

      obs.cv <-
          bru_obs(
            formula = formula(paste0("exposed24 ~ Intercept + Village + hhID + field +", modelSpecsIncidence[which(names(modelSpecsIncidence)==fit.name)])),
            family = "binomial",
            data = dat.all.cv)


      # Fit the model
      cv.fit.start <- Sys.time()
      fit.cv <-
        bru(
          components = components.cv,
          obs = obs.cv,
          options = inla.opt)
      cv.fit.finish <- Sys.time()
      print(cv.fit.finish - cv.fit.start)


      dat.all.inc$pred.fold[which(!cv$biomod_table[, i])] <-
        predict(fit.cv,
              formula = paste0("Intercept + ", modelSpecsIncidence[[fit.name]]),
              newdata = dat.all.cv)$mean[which(!cv$biomod_table[, i])]
      print(table(!is.na(dat.all.inc$pred.fold)))
  }

    #' Back-transform from the cloglog-scale linear predictor to a probability
    dat.all.inc$probPred <- 1 - exp(-exp(dat.all.inc$pred.fold))
    plot(dat.all.inc$probPred, dat.all.inc$exposed24)

    #' exposed24 is NA for animals already seropositive in 2023, so exclude
    #' those from the Brier score (BrierScore() requires strictly binary input)
    scoredRows <- !is.na(dat.all.inc$exposed24)
    brierScore <- BrierScore(dat.all.inc$exposed24[scoredRows], pred = dat.all.inc$probPred[scoredRows])

    #' Brier score per fold, so we can check afterwards whether the overall
    #' ranking is consistent across folds or driven by just one or two
    foldBrier <-
      sapply(1:cluster.k, function(k) {
        rows <- which(dat.all.inc$fold == k & !is.na(dat.all.inc$exposed24))
        BrierScore(dat.all.inc$exposed24[rows], pred = dat.all.inc$probPred[rows])
      })
    cv.diagnostics[[fit.name]] <<- list(probPred = dat.all.inc$probPred, foldBrier = foldBrier)

    brierScore
  })

##### CV diagnostic checks #####

#' 1. Do the models' out-of-fold predicted probabilities actually differ?
#' If they're near-identical, a lower Brier score for the simpler model isn't
#' meaningful - it just means the extra covariates aren't changing anything.
#' i.e. if their correlation is qutie good then the brier score doesn't mean
#' anything.
predMatrix2024 <- sapply(cv.diagnostics, function(x) x$probPred)
cat("Correlation between models' out-of-fold predicted probabilities:\n")
print(cor(predMatrix2024, use = "pairwise.complete.obs"))

#' 2. Is the Brier score ranking consistent across all folds, or is it
#' driven by one or two folds?
#'
foldBrierTab2024 <- sapply(cv.diagnostics, function(x) x$foldBrier)
rownames(foldBrierTab2024) <- paste0("fold", 1:cluster.k)
cat("\nBrier score by fold and model:\n")
print(round(foldBrierTab2024, 3))

cat("\nHow many folds does each model 'win' (lowest Brier)?\n")
foldWinnersInc <- factor(apply(foldBrierTab2024, 1, function(x) names(x)[which.min(x)]),
                      levels = colnames(foldBrierTab2024))
print(table(foldWinnersInc))


#' ## Tabulate results for all models
#' Get WAIC for all models

WAIC.tab.incidence <- cbind(Fixed_effects = unlist(modelSpecsIncidence), IC.inla(mods$outputs[4:length(names(mods$outputs))], crit = "waic"))

covariates <-
   c(
    Culex = "culexAbundance.sds",
    Aedes = "aedesAbundance.sds" )

for(i in 1:nrow(WAIC.tab.incidence)) {
  for(j in 1:length(covariates)) {
    WAIC.tab.incidence[i, "Fixed_effects"] <-
      gsub(covariates[j], names(covariates)[j], WAIC.tab.incidence[i, "Fixed_effects"])
  }
  WAIC.tab.incidence[i, "Fixed_effects"] <- gsub("\\.sds", "", WAIC.tab.incidence[i, "Fixed_effects"])
}
names(WAIC.tab.incidence) <- gsub("delta", "Δ", names(WAIC.tab.incidence))

#' Add 5-fold CV results to WAIC table

WAIC.tab.incidence <- cbind(WAIC.tab.incidence, as.matrix(round(cv.results, 4))[rownames(WAIC.tab.incidence), ])
colnames(WAIC.tab.incidence) <- c("Fixed_effects"
                        ,"WAIC"
                        ,"ΔWAIC"
                        ,"Weights"
                        ,"BrierScore")

write.csv(WAIC.tab.incidence,  "FOI_analysis/public_outputs/modelFitTabIncidence.csv")

##### Extract the posteriors to get Village-level 2024 FOI #####

Xinc <- as.data.frame(model.matrix(~ prop_year +
                                  culexAbundance.sds +
                                    aedesAbundance.sds,
                                data = cattleDatInc))


# extract FOI using the same X object as in the v9 script with aedes and culex and age.
#' Evaluate the spatial field (object field.est) from our model only at our sampled locations
spatialRE24 <-
  fm_evaluate(fm_evaluator(mods$inputs$spatial24.null$mesh, loc = cattleDatInc),
              mods$outputs$spatial24.null$summary.random$field[, point.est])

#' Sum the fixed effects predictions, the village random effect, and the spatial
#' random effect, giving village-level log FOI predictions
fixedPred24 <-
  (as.matrix(Xinc[, 1]) %*% mods$outputs$spatial24.null$summary.fixed[, point.est])
villageRE24 <-
  mods$outputs$spatial24.null$summary.random$Village[
    match(cattleDatInc$Village, mods$outputs$spatial24.null$summary.random$Village$ID), point.est]
cattleDatInc$logvillageFOI24 <-
  fixedPred24 +
  villageRE24 +
  spatialRE24

#' Add the herd random effect, giving herd-level log FOI predictions
herdRE24 <-
  mods$outputs$spatial24.null$summary.random$hhID$`0.5quant`[
    match(as.character(cattleDatInc$hhID), mods$outputs$spatial24.null$summary.random$hhID$ID)]
cattleDatInc$logherdFOI24 <- cattleDatInc$logvillageFOI24 + herdRE24

#' Exponentiate log FOI to FOI
cattleDatInc$villageFOI24 <- exp(cattleDatInc$logvillageFOI24)
cattleDatInc$herdFOI24 <- exp(cattleDatInc$logherdFOI24)

#' This is the RAW foi. What this is basically saying, is an animal would experienced
#' FOI number of of infectious exposure events per year if it kept being exposed indefinitely
#' but this isn't possible because seroconverting is a one time event (in this framewokr)
#' so basically how many exposure events happen on average...

data.frame(
  quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
  villageFOI24 = quantile(cattleDatInc$villageFOI24, probs = c(0, 0.025, 0.5, 0.975, 1)),
  row.names = NULL)

#' But what we want to know actually is what is the chance at least one
#' event happens before the year is out?
#' so this is where this comes in....

#' Annual % of susceptible animals expected to seroconvert (bounded 0-100%,
#' unlike the raw FOI rate above): 1 - exp(-FOI), expressed as a percentage.
#' although it is actually the probability of seroconverting per year.
cattleDatInc$annualPercentSeroconversion <- (1-exp(-cattleDatInc$villageFOI24))*100

data.frame(
  quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
  annualvillageSC = quantile(cattleDatInc$annualPercentSeroconversion, probs = c(0, 0.025, 0.5, 0.975, 1)),
  row.names = NULL)

#' Same summary, broken down by District
cattleDatInc %>%
  st_drop_geometry() %>%
  group_by(District) %>%
  summarise(
    n = n(),
    `0%` = quantile(annualPercentSeroconversion, 0),
    `2.5%` = quantile(annualPercentSeroconversion, 0.025),
    `50%` = quantile(annualPercentSeroconversion, 0.5),
    `97.5%` = quantile(annualPercentSeroconversion, 0.975),
    `100%` = quantile(annualPercentSeroconversion, 1)
  )

#' at low FOIs FOI ~= to the annual percent, but they diverge at higher FOIs.
#' so an expected count (i.e. incident rate)

#' Express FOI as a standard Poisson incidence rate (cases per 100 animal-years
#' at risk). FOI is already a per-animal-year rate because exposed.yrs was
#' used as the model offset, so this is a re-expression, not a new estimate:
#' N infections in T animal-years of exposure ~ Poisson(FOI * T)
incidenceRate24.tab <-
  data.frame(
    quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
    incidenceRate100AY = quantile(cattleDatInc$villageFOI24, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = TRUE) * 100,
    row.names = NULL)

write.csv(incidenceRate24.tab, "FOI_analysis/public_outputs/incidenceRatePer100AnimalYearsIncidence.csv", row.names = FALSE)

VillagePred$FOI24 <-
  tapply(cattleDatInc$villageFOI24, cattleDatInc$Village, mean)[as.character(VillagePred$Village)]

VillagePred$annualPercentSeroconversion2024 <- (1 - exp(-VillagePred$FOI24)) * 100
VillagePred$timesHigher <-
  round(VillagePred$annualPercentSeroconversion2024/VillagePred$annualPercentSeroconversion)

write.csv(VillagePred, "FOI_analysis/public_outputs/villagePred.csv", row.names = F)

##### Main figure: how many times higher #####

(changePlot <-
  ggplot() +
  annotation_map_tile(type = "cartolight", zoomin = 0) +  # You can change the tile type
  geom_sf(data = VillagePred,
          aes(fill = timesHigher,
              shape = District),
          colour = "darkgrey",
          size = 4) +
  scale_fill_viridis_c(option = "D",
                       direction = -1,
                       breaks = scales::breaks_width(2)) +
  scale_shape_manual(values = c(Babati = 21, Hai = 22, Moshi = 24)) +
  annotation_scale(location = "tl", width_hint = 0.5) +
  coord_sf(crs = crs.km) +
  theme_minimal() +
   theme(
     legend.key.size = unit(0.8, "cm"),
     legend.spacing.y = unit(0.3, "cm"),
     plot.margin = margin(2, 2, 2, 2),
     legend.position = "bottom",
     legend.title = element_text(size = 12)
   )+
  labs(fill = "x-fold increase in \nannual seroconversion \n2023-2024")
  )

cowplot::plot_grid(changePlot,
          vectorIncidence,
          nrow = 1, ncol = 2, labels = c("A", "B")
)
ggsave("FOI_analysis/public_outputs/FOIpanelIncidence.pdf")
