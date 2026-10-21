# Public version of genus-level-spatial-analysis.R, run on the anonymised data in data_public/.

# ## BACKGROUND ----

# The aim of this script is to fit a spatial GLMM to the 
# NERC RVF project vector trapping data,
# pooling vector species at the genus level, and use it to predict vector abundance.
# This will be possible if we find 
# predictors / covariates that strongly predict vector abundance. 
# The vector abundance predictions will then be used to predict 
# cattle force of infection (FOI).
#   
# There is trapping data on three genera, Culex, Aedes and Mansonia, 
# but only Culex and Aedes were trapped in sufficient numbers.  
 
# This script gratefully used or adapted code from 
# Greg Albery's INLA tutorial:
# https://ourcodingclub.github.io/tutorials/inla/
# and uses some of the plotting and summary functions in his package.  

# Other useful tutorials:  
# https://datascienceplus.com/spatial-regression-in-r-part-2-inla/  
# https://becarioprecario.bitbucket.io/inla-gitbook/  
#   
#   

# ## PACKAGES AND FUNCTIONS ----

# Load packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(magrittr)
  library(tidyverse)
  library(RColorBrewer)
  library(beepr)
  library(ggmap)
  library(ggspatial)
  library(ggtext)
  library(fmesher)
  library(scales)
  library(sf)
  library(MASS) # for glm.nb, ML GLM for comparison with INLA
  library(lme4) # for ICC calculation
  library(DHARMa)
  library(MCMCglmm); library(coda); library(ape) # required by the Efxplot function
  library(prettymapr)
  library(parallel)
  library(knitr)
  library(patchwork)
  library(colorspace)
  library(formattable)
  library(htmltools)
  library(htmlwidgets)
  library(blockCV)
  library(RANN) # for nn2, nearest neighbour function
  #library(ggregplot) # devtools::install_github("gfalbery/ggregplot")
  # To install R-INLA for the first time:
  # install.packages("INLA",
  #                  repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"),
  #                  dep=TRUE)
  # Check for upgrades:
  # INLA::inla.upgrade() 
  # NOTE on 25 Aug 2026: fitting the models became impossibly slow after 
  # updating to the current version of INLA (26.08.07), which was resolved after 
  # reverting to version 25.10.19.
  library(INLA)
  library(inlabru)
})

# Clear the global environment and close graphics devices
rm(list = ls())
graphics.off()

# Options
options(warn = 1)

# Start timer
global.start <- Sys.time()


# Set random seed for reproducibility. All random seeds will be chosen using
# https://www.random.org/integers/?num=1&min=0&max=1000000000&col=1&base=10&format=html&rnd=new
set.seed(790092629)

# # Set INLA to use one core to prevent interference with mclapply
# inla.setOption(num.threads = "1:1", inla.mode = "classic")

# Load local functions
source("functions.R")

dir.create("figures_public", showWarnings = FALSE)
dir.create("Vector_analysis/public_outputs/fixed_effects_estimates", recursive = TRUE, showWarnings = FALSE)

# ## GLOBAL ANALYSIS CHOICES ----

# Which vector genus to analyse? The following genera are available:
genera <- c("Culex", "Aedes", "Mansonia")

# Choose colours for the three genera (from the Okabe–Ito palette):
genera.cols <- c("#E69F00", "#56B4E9", "#009E73")
names(genera.cols) <- genera

# Only Culex and Aedes have enough numbers for modelling abundance. 
genus <- genera[2] # Choose genus (1: Culex or 2: Aedes)

# ## LOAD AND PROCESS DATA ----

# ### Load input datasets

# Load the genus-level vector abundance data. This data set was created by the script 
# "./Dropbox/RVFV_data/mosquitoes/dataAndRCode/r1_readAndCleanData.R" and written to CSV using:
# write.csv(mos, "GenusLevel.csv", row.names = FALSE)
mos.raw <- read.csv("data_public/GenusLevel.csv")
dim(mos.raw)
names(mos.raw)

# Load the data set with trap-level covariates
trap <- read.csv("data_public/trapSetWithCovariatesFeb2025.csv")
dim(trap)
names(trap)

# Load the rainfall data
rain <- read.csv("data_public/rainfallDataSummaries.csv")
rain$X <- NULL
dim(rain)
names(rain)

# Load the cattle serology data and the covariates, including lifetime
# NDVI and rain days summaries for cattle sampling sites
cattle <- read.csv("data_public/allCattleDat.csv")
dim(cattle)

# Drop any animals that have no serology result from 2023
cattle <- cattle[!is.na(cattle$status2023), ]
dim(cattle)

# Load the cattle incidence data
cattle.inc <- read.csv("data_public/cattleIncidenceData.csv")
dim(cattle.inc)
names(cattle.inc)

# List of all quarters covered by the lifetimes of the cattle
all.quarters <- expand.grid(quarter = paste0("Q", 1:4), year = 2023:2008, stringsAsFactors = FALSE)
all.quarters$yq <- paste(all.quarters$year, all.quarters$quarter, sep = "_")

# Missing NDVI and rain were imputed before anonymisation
mean(is.na(cattle[, c("lifetimeNdvi500m.mean", "lifetimeNdvi500m.sd", "lifetimeLog10RainDays500m.mean")]))
mean(is.na(cattle.inc[, c("mean_ndvi_500m", "mean_rain_days_500m")]))

# ### Harmonise identifiers and variable names

# #### Harmonise trap identifiers

# Paste "TS" onto trap ID number to avoid any numeric-vs-character indexing errors
rain$TS_ID <- paste0("TS", rain$TS_ID)
trap$TS_ID <- paste0("TS", trap$TS_ID)

# Does each trap have a single row of covariate values?
length(trap$TS_ID) == length(unique(na.omit(trap$TS_ID)))
# ...yes. We can use TS_ID as row names:
rownames(trap) <- trap$TS_ID

# #### Harmonise cattle variable names

# Rename variables in the data frame 'cattle' to match the other data sets
cattle$elevation500m <- cattle$elev500m
cattle$popDens500m <- cattle$pd500m
cattle$buildingDensity500m <- cattle$bd500m
cattle$proportionIrrigatedCrop <- cattle$propI500m
cattle$elev500m <- cattle$pd500m <- cattle$bd500m <- cattle$propI500m <- NULL

# Create cattle birth year variable
cattle$birthyear <- 2023 - cattle$newAgeInt

# Village names are already coded in data_public/

# #### Harmonise factor coding

# Remove the hyphen from peri-urban (it interferes with fitting 
# Class via a model matrix)
trap$Class <- gsub("-", "", trap$Class)

# #### Add the rain data to the trap covariate data

# Which variables are in both trap covariate and rain data sets?
intersect(names(trap), names(rain))

# Which variables are unique to the trap covariate data set?
setdiff(names(trap), names(rain))

# Which variables are unique to the rain data set?
setdiff(names(rain), names(trap))

# Add the rain variables to the trap data set, then 
# remove the rain data 
trap$meanRainfall <- rain$meanRainfall[match(trap$TS_ID, rain$TS_ID)]
trap$NRainDays <- rain$NumberRainDays[match(trap$TS_ID, rain$TS_ID)]
rm(rain)

# Which variables are in both trap covariate and vector data sets?
intersect(names(trap), names(mos.raw))

# Which variables are unique to the covariate data set?
setdiff(names(trap), names(mos.raw))

# Which variables are unique to the vector data set?
setdiff(names(mos.raw), names(trap))

# ### Merge mosquito abundance and covariate datasets

# #### Standardise mosquito identifiers

# Paste TS onto trap ID number to avoid any numeric/character indexing errors
mos.raw$TS_ID <- paste0("TS", mos.raw$TS_ID)

# There are two columns named "...1" and "...10"
# Remove them:
mos.raw$...10 <- mos.raw$...1 <- NULL

# Which variables are in both trap covariate and vector data sets?
intersect(names(trap), names(mos.raw))

# Which variables are unique to the covariate data set?
setdiff(names(trap), names(mos.raw))

# Which variables are unique to the vector data set?
setdiff(names(mos.raw), names(trap))

# We should take all the covariate data from the covariate data set, so 
# remove the covariate columns in the vector data, keeping only the unique columns,
# and TS_ID which we will use to link to the covariate data
mos <- mos.raw[, c("TS_ID", setdiff(names(mos.raw), names(trap)))]
dim(mos)

# Check there is one row per trap in the vector trapping data
length(unique(na.omit(mos$TS_ID))) == nrow(mos)

# There is, so we can use trap ID as rownames
rownames(mos) <- mos$TS_ID

# Merge the covariates onto the vector data. First check that all the traps 
# in the vector abundance data set are in the covariates data set:
setdiff(mos$TS_ID, trap$TS_ID)

# "character(0)" means they are all present. Are there traps with covariate data
# but no vector trapping data? (Not a problem if there are.)
setdiff(trap$TS_ID, mos$TS_ID)

# First change the name of the covariate data trap ID, so we don't have columns
# with the same name
trap$TS_ID.trap <- trap$TS_ID
trap$TS_ID <- NULL

# Merge the two data frames
mos <- cbind(mos, trap[as.character(mos$TS_ID), ])
sort(names(mos))

# Check that the trap IDs match up
all(mos$TS_ID == mos$TS_ID.trap)
mos$TS_ID.trap <- NULL

# #### Create analysis factors

# Now each trap count has linked covariate data.   

# District names.
# There is only one village in Moshi district, which will look odd on some plots,
# so create a 4-district factor that combines Moshi with neighbouring Hai.
mos$District4 <- 
  factor(mos$District, 
         c("Babati", "Hai", "Monduli", "Moshi", "Same"), 
         c("Babati", "Hai/Moshi", "Monduli", "Hai/Moshi", "Same"))
table(mos$District4, exclude = NULL)

# Create a 5-district factor
mos$District <- factor(mos$District)
table(mos$District, exclude = NULL)

# Sort mos by district then village
mos <- mos[order(mos$District, mos$Village), ]

# Make factors for district, village, trap, sampling round, class (peri-urban vs rural)

# Sampling round
mos$samplingRound <- factor(gsub("-", "", (mos$samplingRound)), c("Mar23", "Aug23", "May24"))

# Peri-urbanicity class
mos$Class <- factor(mos$Class, c("Rural", "Periurban"))

# District and village
mos$District <- factor(mos$District)
table(mos$District, exclude = NULL)
mos$Village <- factor(mos$Village)
table(mos$Village, exclude = NULL)

# Some villages were sampled once (Cross-sectional), others over three rounds.  
# Create a factor to record this
VillageTimesSampled <- rowSums(table(mos$Village, mos$samplingRound) > 0)
mos$TimesSampled <- factor(VillageTimesSampled[as.character(mos$Village)])
table(mos$TimesSampled, exclude = NULL)

# Time-varying village random effect factor, levels grouped by sampling round
mos$VillageTV <- 
  factor(paste(mos$Village, "-", mos$samplingRound),
         apply(expand.grid(levels(mos$Village), levels(mos$samplingRound)), 1, paste, collapse = " - "))
mos$VillageTV <- droplevels(mos$VillageTV)
table(mos$VillageTV, exclude = NULL)

# Time-varying village random effect factor, levels grouped by village
mos$VillageTV.round.order <- 
  factor(paste(mos$samplingRound, "-", mos$Village),
         apply(expand.grid(levels(mos$samplingRound), levels(mos$Village)), 1, paste, collapse = " - "))
mos$VillageTV.round.order <- droplevels(mos$VillageTV.round.order)
table(mos$VillageTV.round.order, exclude = NULL)

# ### Create mosquito analysis variables

# #### First explore relationships among genera

# Are Culex and Aedes correlated?
cor(mos$Culex, mos$Aedes, method = "spearman")
plot(1 + mos$Culex, 1 + mos$Aedes, log = "xy")
cor(1 + tapply(mos$Culex, mos$Village, mean), 1+ tapply(mos$Aedes, mos$Village, mean), 
    method = "spearman")
plot(1 + tapply(mos$Culex, mos$Village, mean), 1+ tapply(mos$Aedes, mos$Village, mean),
     log = "xy")

# #### Format identifiers and dates

# Make trap ID a factor
mos$TS_ID <- factor(mos$TS_ID)

# Convert trapping dates to date format
mos$Date_deployed <- as.Date(mos$Date_deployed, "%d/%m/%Y")
mos$Date_collected <- as.Date(mos$Date_collected, "%d/%m/%Y")

# For how many days were traps deployed?
table(mos$Date_collected - mos$Date_deployed)

# ### Evaluate trap types and justify exclusions

# #### Create trap-type factor

# Type of trap used
mos$Trap_type <- factor(mos$Trap_type, c("CDC Light Trap", "BG Sentinel"), c("CDCLT", "BGS"))
table(mos$Trap_type, exclude = NULL)
table(mos$samplingRound, mos$Trap_type, exclude = NULL)

# 72 of 551 traps deployed were BG Sentinel Traps, the remainder were CDC Light Traps.
# The BGS traps were deployed quite widely in space and over time
table(mos$Trap_type)
table(mos$TimesSampled, mos$Trap_type)
table(mos$District, mos$Trap_type)

# #### Compare vector catches by trap type

# Did the BGS traps catch many vectors?
tapply(mos$Culex, mos$Trap_type, sum)
tapply(mos$Aedes, mos$Trap_type, sum)
tapply(mos$Mansonia, mos$Trap_type, sum)
tapply(mos$Culex, mos$Trap_type, mean)
tapply(mos$Aedes, mos$Trap_type, mean)
tapply(mos$Mansonia, mos$Trap_type, mean)

# What proportion of vectors were trapped by each trap type?
prop.table(tapply(mos$Culex, mos$Trap_type, sum) + tapply(mos$Aedes, mos$Trap_type, sum) + 
             tapply(mos$Mansonia, mos$Trap_type, sum))

# BG sentinels trapped 2 Culex per trap (total 160, 0.6% of the total), 
# and a total of 1 Aedes and 0 Mansonia.   
# Are the trap counts correlated at the village level for Culex?
cor(tapply(mos$Culex, list(mos$Village, mos$Trap_type), mean), 
    use = "pairwise", method = "spearman")
# ...no.  
# Taking together the very low numbers of mosquitoes trapped by BGS and the lack of correlation
# with CDCLT, it seems very unlikely that BGS will contribute so let's drop these traps.
mos <- droplevels(mos[mos$Trap_type %in% "CDCLT", ])
nrow(mos)

# #### Exclude BG Sentinel traps and summarise sampling

# Duration of each sampling round
sampling.dates.tab <-
  cbind(
    `Sampling round` = levels(mos$samplingRound),
    do.call("rbind", tapply(mos$Date_deployed, mos$samplingRound, function(x) as.character(range(x)))),
    `Duration (nights)` = tapply(mos$Date_deployed, mos$samplingRound, function(x) diff(range(x))) + 1,
    `N villages` = tapply(mos$Village, mos$samplingRound, function(x) length(unique(x))))
colnames(sampling.dates.tab)[2:3] <- c("First deployment date", "Last deployment date")
rownames(sampling.dates.tab) <- NULL
formattable(as.data.frame(sampling.dates.tab))
# ...approx 1 month per round.  

# Number of traps by district and sampling round
rbind(cbind(table(mos$District, mos$samplingRound), Total = table(mos$District)),
      Total = c(table(mos$samplingRound), nrow(mos)))

# Numbers trapped of each genus:
colSums(mos[, genera])

# Total number trapped:
sum(colSums(mos[, genera]))

# Mean numbers trapped per sampling round
tapply(mos$Culex, mos$samplingRound, mean)
tapply(mos$Culex, mos$samplingRound, sd)
tapply(mos$Aedes, mos$samplingRound, mean)
tapply(mos$Aedes, mos$samplingRound, sd)
tapply(mos$Culex > 0, mos$samplingRound, mean)
tapply(mos$Aedes > 0, mos$samplingRound, mean)
mean(mos$Culex > 0)
mean(mos$Aedes > 0)

# Store trap counts of the chosen genus as Total.Females, which will
# be the response variable in the GLMMs
mos$Total.Females <- mos[, genus]

# Mean number of females per trap
mean(mos$Total.Females)

# ### Create environmental and demographic covariates

# #### Create buidling and population density categories

# Create a factor that divides locations into high/low building density
# combined with high/low population density
mos$buildingDensity500m.cat <- 
  cut(mos$buildingDensity500m, c(-0.0001, median(mos$buildingDensity500m), Inf),
      labels = c("Low", "High"))
mos$popDens500m.cat <- 
  cut(mos$popDens500m, c(-0.0001, median(mos$popDens500m), Inf),
      labels = c("Low", "High"))
table(mos$buildingDensity500m.cat, mos$popDens500m.cat)

# Factor combining these two dichotomised factors
mos$popDensbuildDens500m.cat <- 
  factor(paste0(mos$buildingDensity500m.cat, mos$popDens500m.cat),
         c("LowLow", "LowHigh", "HighLow", "HighHigh"),
         c("LoBD-LoPD", "LoBD-HiPD", "HiBD-LoPD", "HiBD-HiPD"))
table(mos$popDensbuildDens500m.cat)


# #### Explore density relationships

# Building density and population density are correlated, although the nature 
# of the correlation depends on district
ggplot(mos, aes(y = buildingDensity500m, x = popDens500m, 
                shape = popDensbuildDens500m.cat, color = District)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = coef(lm(buildingDensity500m ~ -1 + popDens500m, data = mos)),
              intercept = 0)
cor(mos$popDens500m, mos$buildingDensity500m)


# ...therefore create a population density per building density measure
mos$popDens500m.per.bd <- mos$popDens500m/mos$buildingDensity500m
ggplot(mos, aes(x = buildingDensity500m, y = popDens500m.per.bd, shape = District, color = District)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", se = FALSE)
cor(mos$popDens500m.per.bd, mos$buildingDensity500m)
# Now the variables are negatively correlated, although the correlation is slightly weaker.  
#   

# ### Select and transform covariates

# #### Define candidate predictors

# Which covariates do we think might predict vector abundance?
all.covariates <- 
  c("samplingRound", 
    "ndvi", 
    "NRainDays", 
    "elevation500m",
    "buildingDensity500m", 
    "popDens500m", 
    "popDens500m.per.bd",
    "popDensbuildDens500m.cat",
    "proportionIrrigatedCrop",
    "proportionCropLandC")

# #### Assess variable distributions

# Check for skew in continuous variables
cont.cov <- all.covariates[!sapply(mos[, all.covariates], is.factor)]
par(mfrow = c(ceiling(sqrt(length(cont.cov))), ceiling(sqrt(length(cont.cov)))))
invisible(sapply(cont.cov, function(cc) hist(mos[, cc], main = cc)))
par(mfrow = c(ceiling(sqrt(length(cont.cov))), ceiling(sqrt(length(cont.cov)))))
invisible(sapply(cont.cov, function(cc) hist(log10(mos[, cc]), main = paste("log10", cc))))
par(mfrow = c(1, 1))

# log10 transform the following variables, to reduce skew, 
# drawing in potentially influential data points:
skewed.vars <-
  c("NRainDays", "buildingDensity500m", "popDens500m", "popDens500m.per.bd")
for(sv in skewed.vars) {
  mos[, paste0(sv, ".log10")] <- log10(mos[, sv])
}


# #### Define modelling covariates

# Re-select the covariates choosing the logged versions
covariates <- 
  c(NDVI = "ndvi", 
    `n rainy days` = "NRainDays.log10", 
    `proportion irrigated` = "proportionIrrigatedCrop",
    elevation = "elevation500m",
    `building density` = "buildingDensity500m.log10", 
    `population density` = "popDens500m.log10") 

# How much variation within and between villages in the continuous covariate values?  
# Graphically and using ICC:
cont.cov <- covariates[!sapply(mos[, covariates], is.factor)]
if(length(cont.cov) > 0) {
  old.par <- 
    par(mfrow = c(ceiling(sqrt(length(cont.cov))), floor(sqrt(length(cont.cov)))),
        tcl = -0.2,
        mgp = c(3, 0.3, 0),
        mar= c(6.3, 4.1, 3.6, 2.1))
  sapply(cont.cov, 
         function(cc) {
           plot(formula(paste(cc, "~ VillageTV.round.order")), data = mos, las = 3,
                xlab = "", cex.axis = 0.6)
           fit <- lmer(paste(cc, "~ (1 | samplingRound + VillageTV + Village)"), data = mos)
           print(cc)
           print(unlist(VarCorr(fit)))
           vill.var <- VarCorr(fit)$Village
           villtv.var <- VarCorr(fit)$VillageTV
           round.var <- VarCorr(fit)$samplingRound
           resid.var <- attr(VarCorr(fit), "sc")^2
           tot.var <- vill.var + villtv.var + round.var + resid.var
           icc.tot <- (vill.var + villtv.var + round.var) / tot.var
           icc.vill <- vill.var / tot.var
           title(paste0(cc, "\nInter-village ICC = ", round(100 * icc.vill), 
                        "% (total ICC = ", round(100 * icc.tot),
                        "%)"))
         })
  par(old.par)
  
  # How closely correlated are the covariates
  panel.hist <- function(x, ...) {
    usr <- par("usr")
    old.par <- par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, ...)
    par(old.par)
  }
  
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    old.par <- par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = 2)
    par(old.par)
  }
  pairs(mos[, cont.cov], upper.panel = panel.cor, diag.panel = panel.hist,
        gap=0, row1attop=FALSE, cex.labels = 1.3)
  
  # No excessively large correlations
}

# #### Standardise continuous predictors

# Convert all continuous covariates to standard deviation scores, mainly to 
# centre and standardise the priors on the intercept and slopes
for(cc in covariates) {
  if(!is.factor(mos[, cc])) {
    mos[, paste0(cc, ".sds")] <- scale(mos[, cc])
    covariates[covariates == cc] <- paste0(cc, ".sds")
  }
}; rm(cc)

# ### Create cattle covariates on model scale

# #### Scale incidence predictors

# In the cattle incidence data, create SD score versions of mean NDVI, NDVI^2 and 
# no of rain days 
cattle.inc$ndvi.sds <- 
  (cattle.inc$mean_ndvi_500m - attr(mos$ndvi.sds, "scaled:center")) / 
  attr(mos$ndvi.sds, "scaled:scale")
cattle.inc$ndvi.sds2 <- cattle.inc$ndvi.sds^2
cattle.inc$NRainDays.log10.sds <- 
  (log10(cattle.inc$mean_rain_days_500m) - attr(mos$NRainDays.log10.sds, "scaled:center")) / 
  attr(mos$NRainDays.log10.sds, "scaled:scale")
cattle.inc$NRainDays.log10.sds2 <- cattle.inc$NRainDays.log10.sds^2

# Add other covariates to the cattle incidence data
for(v in c("elevation500m", "buildingDensity500m.log10", 
           "popDens500m.log10", "proportionIrrigatedCrop")) {
  if(substr(v, nchar(v) - 5, nchar(v)) == ".log10") {
    cattle.inc[, paste0(v, ".sds")] <- 
      (log10(cattle[match(cattle.inc$newID, cattle$newID), substr(v, 1, nchar(v) - 6)]) - attr(mos[[paste0(v, ".sds")]], "scaled:center")) / 
      attr(mos[[paste0(v, ".sds")]], "scaled:scale")
  } else {
    cattle.inc[, paste0(v, ".sds")] <- 
      (cattle[match(cattle.inc$newID, cattle$newID), v] - attr(mos[[paste0(v, ".sds")]], "scaled:center")) / 
      attr(mos[[paste0(v, ".sds")]], "scaled:scale")
  }
}

# Make cattle incidence data variable for the interaction between irrigation and rain 
cattle.inc$IrrigatedXRain <- cattle.inc$proportionIrrigatedCrop.sds * cattle.inc$NRainDays.log10.sds

# Create NDVI, NDVI^2 and no of rain days variables averaged across the
# lifespan of each animal, from the lifetime mean and SD of NDVI and
# lifetime mean log10 rain days
ndvi.c <- attr(mos$ndvi.sds, "scaled:center")
ndvi.s <- attr(mos$ndvi.sds, "scaled:scale")
cattle$ndvi.sds <- (cattle$lifetimeNdvi500m.mean - ndvi.c) / ndvi.s
cattle$ndvi.sds2 <-
  ((cattle$lifetimeNdvi500m.mean - ndvi.c)^2 + cattle$lifetimeNdvi500m.sd^2) / ndvi.s^2
cattle$NRainDays.log10.sds <-
  (cattle$lifetimeLog10RainDays500m.mean - attr(mos$NRainDays.log10.sds, "scaled:center")) /
  attr(mos$NRainDays.log10.sds, "scaled:scale")

# #### Create cattle location summaries

# Create village level coordinates for cattle sampling locations
cattle$VLongitude <- tapply(cattle$long, cattle$Village, mean)[cattle$Village]
cattle$VLatitude <- tapply(cattle$lat, cattle$Village, mean)[cattle$Village]

# Make peri-urban classification match mosquito data set
cattle$classification <- gsub("-", "", cattle$classification)

# District and village factors
cattle$District <- factor(cattle$District)
cattle$Village <-  factor(cattle$Village)

# Create empty (NA) fields for rain and NDVI data in the cattle data frame
cattle$NRainDays5y.log10.sds <- NA
cattle$NRainDays8m.log10.sds <- NA
cattle$ndvi5y.dry.sds <- NA
cattle$ndvi5y.sds <- NA
cattle$ndvi8m.sds <- NA


# To what extent do villages overlap between mos and cattle?
setdiff(mos$Village, cattle$Village)
setdiff(cattle$Village, mos$Village)
intersect(cattle$Village, mos$Village)
setdiff(mos$District, cattle$District)
setdiff(cattle$District, mos$District)
intersect(cattle$District, mos$District)

# Sort cattle by district then village
cattle <- cattle[order(cattle$District, cattle$Village), ]

# Impute missing coordinates as the mean of the non-missing coordinates in
# the same village
sum(is.na(cattle[, c("lat", "long")]))
table(rowSums(is.na(cattle[, c("lat", "long")])), cattle$Village)
for(vill in levels(cattle$Village)) {
  cowvill <- cattle[cattle$Village == vill, ]
  cowvill$lat[is.na(cowvill$lat)] <- mean(cowvill$lat, na.rm = TRUE)
  cowvill$long[is.na(cowvill$long)] <- mean(cowvill$long, na.rm = TRUE)
  cattle[cattle$Village == vill, c("lat", "long")] <- cowvill[, c("lat", "long")]
  rm(cowvill)
}
table(rowSums(is.na(cattle[, c("lat", "long")])), cattle$Village)


# Impute missing building density as the mean of the rest of the village
if(length(unique(cattle$Village[is.na(cattle$buildingDensity500m)])) == 1) {
  na.builddens.village <- unique(cattle$Village[is.na(cattle$buildingDensity500m)])
  cattle$buildingDensity500m[is.na(cattle$buildingDensity500m)] <-
    mean(cattle$buildingDensity500m[cattle$Village %in% na.builddens.village], na.rm = TRUE)
}

# #### Finalise cattle outcome dataset

# Restrict the cattle serology data set to rows with serology data from 2023
cattle <- droplevels(cattle[!is.na(cattle$status2023), ])
dim(cattle)

# Classify positives as exposed and doubtfuls and negatives as not exposed
cattle$exposed <- as.integer(cattle$status2023 == "Positive")
table(cattle$exposed, exclude = NULL)

# ### Merge the vector and livestock data sets

# #### Create cattle prediction dataset

# Create a mosquito vector prediction data set, mosp, which will store the coordinates
# and fixed and random effect data to allow prediction of vector abundance to locations
# where cattle serology was recorded.
mosp <- mos[0, ]
mosp[1:nrow(cattle), ] <- NA
dim(mosp)

# Create columns in mosp to match mos
mos$newID <- NA
mosp$newID <- cattle$newID
mos$hhID <- NA
mosp$hhID <- cattle$hhID
mos$newAgeInt <- NA
mosp$newAgeInt <- cattle$newAgeInt
mos$birthyear <- NA
mosp$birthyear <- cattle$birthyear
mos$exposed <- NA
mosp$exposed <- cattle$exposed
mosp$Longitude <- cattle$long
mosp$Latitude <- cattle$lat
mosp$VLongitude <- cattle$VLongitude
mosp$VLatitude <- cattle$VLatitude
mosp$Village <- cattle$Village
mosp$District <- cattle$District
mosp$samplingRound[1:nrow(mosp)] <- "Mar23"
mosp$TimesSampled[1:nrow(mosp)] <- "1"
mosp$VillageTV <- 
  factor(paste(mosp$Village, "-", mosp$samplingRound),
         apply(expand.grid(levels(mosp$Village), levels(mosp$samplingRound)), 1, paste, collapse = " - "))
mosp$VillageTV <- droplevels(mosp$VillageTV)
mosp$Trap_type[1:nrow(mosp)] <- "CDCLT"
mosp$Class[1:nrow(mosp)] <- cattle$classification
mosp$elevation500m.sds <- 
  (cattle$elevation500m - attr(mos$elevation500m.sds, "scaled:center")) / 
  attr(mos$elevation500m.sds, "scaled:scale")
mosp$buildingDensity500m.log10.sds <- 
  (log10(cattle$buildingDensity500m) - attr(mos$buildingDensity500m.log10.sds, "scaled:center")) / 
  attr(mos$buildingDensity500m.log10.sds, "scaled:scale")
mosp$popDens500m.log10.sds <- 
  (log10(cattle$popDens500m) - attr(mos$popDens500m.log10.sds, "scaled:center")) / 
  attr(mos$popDens500m.log10.sds, "scaled:scale")
mos$popDens500m.per.bd.log10.sds <- scale(log10(mos$popDens500m/mos$buildingDensity500m))
mosp$popDens500m.per.bd.log10.sds <- 
  (log10(cattle$popDens500m/cattle$buildingDensity500m) - 
     attr(mos$popDens500m.per.bd.log10.sds, "scaled:center")) /
  attr(mos$popDens500m.per.bd.log10.sds, "scaled:scale")
mosp$proportionIrrigatedCrop.sds <- 
  (cattle$proportionIrrigatedCrop - attr(mos$proportionIrrigatedCrop.sds, "scaled:center")) / 
  attr(mos$proportionIrrigatedCrop.sds, "scaled:scale")

# Note that the dynamic variables (NDVI and rain) will be NA or 999 because of the change 
# in the method of predicting vector abundance
mosp$ndvi.sds <- cattle$ndvi.sds 
mosp$NRainDays.log10.sds <- cattle$NRainDays.log10.sds
mosp$ndvi8m.sds <- cattle$ndvi8m.sds # use the 8 month mean for incidence
mosp$NRainDays8m.log10.sds <- cattle$NRainDays8m.log10.sds # use the 8 month mean for incidence

# ### Convert coordinates and create spatial objects

# #### Project coordinates

# Project longitude and latitude onto a plane using the Mercator projection, which gives
# coordinates in metres (locally) by default. I want the units to be km because INLA works
# better with smaller scale numbers. Define a custom CRS using a PROJ string:
crs.km <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
mos <- st_transform(st_as_sf(mos, coords = c("Longitude", "Latitude"), crs = 4326), crs = crs.km)
mosp <- st_transform(st_as_sf(mosp, coords = c("Longitude", "Latitude"), crs = 4326), crs = crs.km)

# Round the eastings and northings to the nearest 10m, as this is close to the limit of GPS accuracy
mos <-
  st_set_geometry(mos,
                  st_as_sf(data.frame(round(st_coordinates(mos), 1)),
                           coords = c("X", "Y"), crs = st_crs(mos))$geometry)
mosp <-
  st_set_geometry(mosp,
                  st_as_sf(data.frame(round(st_coordinates(mosp), 1)),
                           coords = c("X", "Y"), crs = st_crs(mosp))$geometry)

# #### Assess spatial sampling structure

# Note that some trap locations are re-used (this is because I rounded to the nearest 10m)
table(table(apply(st_coordinates(mos$geometry), 1, paste, collapse = "-")))
# How many unique locations are there?
length(unique(apply(st_coordinates(mos$geometry), 1, paste, collapse = "-")))


# How far apart are village centroids?
# Inter-village distances:
villages <- mos %>% 
  group_by(Village) %>% 
  summarise(geometry = st_centroid(st_union(geometry)))
intervillage.km <- as.numeric(dist(st_coordinates(villages)))
summary(intervillage.km)

# How far apart are traps placed in villages?
# Inter-trap distances within villages at each timepoint:
intertrap.km.by.VillageTV <-
  as.numeric(unlist(sapply(levels(mos$VillageTV), function(v) c(dist(st_coordinates(mos[mos$VillageTV == v, ]))))))
summary(intertrap.km.by.VillageTV)

# Histograms of inter-village and inter-trap distances
wrap_plots(
  ggplot() + 
    geom_histogram(aes(x = intervillage.km), fill = "grey",
                   colour = "black") +
    ylab("Frequency") +
    xlab("Inter-village distance (km)") +
    theme_minimal(), 
  ggplot() + 
    geom_histogram(aes(x = intertrap.km.by.VillageTV), fill = "grey",
                   colour = "black") +
    ylab("Frequency") +
    xlab("Inter-trap distance within villages (km)") +
    theme_minimal(),
  ncol = 1)
ggsave("figures_public/intervillage_km.pdf", height = 4, width = 4)

# ### Check for missingness

colSums(is.na(mos))
# Livestock has one missing value, but we are not using livestock data, so 
# no need to impute

# ### Descriptive analyses and exploratory plots

# #### Summarise village abundance

# Some descriptive analysis first. Make a table of mean number of
# adult females trapped, by village and sampling round:
mean.total.by.vill <- 
  tapply(mos$Total.Females, st_drop_geometry(mos)[, c("Village", "samplingRound")], 
         mean)
vill.tab <- table(mos$Village, mos$samplingRound)
vill.tab.final <- 
  matrix(paste0(round(mean.total.by.vill, 1), " (n=", vill.tab, ")"), 
         nrow = nrow(vill.tab), 
         dimnames = dimnames(vill.tab))
vill.tab.final[vill.tab.final == "NA (n=0)"] <- "-"

# The table below shows mean (n=) number of adult females trapped, 
# by village and sampling round:
print(vill.tab.final, quote = FALSE)
formattable(as.data.frame(vill.tab.final))

# ...looking at the villages that are sampled three times, is there a
# consistent village effect? Are there clear differences in abundance 
# between time points? No obvious consistency over time within villages,
# except Babati-Kisangaji at the 1st and 3rd sampling points.  

ggplot(data = mos, mapping = aes(x = Village, y = log10(1 + Total.Females), colour = District)) +
  geom_boxplot() + 
  facet_wrap(~ samplingRound, ncol = 1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# #### Visualise temporal abundance trends

# Make plot showing change in abundance over time 
y.labels <- c(0, 10^(1:4))
abund.over.time <- 
  ggplot(data = mos[mos$TimesSampled == "3", ], 
         mapping = aes(x = samplingRound, y = log10(1 + Total.Females), 
                       colour = Village, group = Village)) +
    geom_jitter(height = 0, width = 0.1, alpha = 0.4) + 
    stat_summary(fun = mean, geom = "point", size = 3, shape = 5) +
    stat_summary(fun = mean, geom = "line") +
    annotate("rect", xmin = 2.2, xmax = 2.8, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "red") +
    annotate("text", x = 2.5, y = Inf, label = "El Niño", vjust = 2, hjust = 0.5, colour = "red") +
    facet_wrap(~ District4, ncol = 1) +
    ylab("N females trapped") +
    scale_y_continuous(
      breaks = log10(1 + y.labels),  # log10(1 + count) values
      labels = y.labels,
      limits = c(0, log10(1 + max(y.labels)))) +
    theme_minimal() +
    ggtitle(paste0("*", genus, "* abundance by district and<br>sampling round (N=", sum(mos$Total.Females), ")")) +
    theme(axis.title.x=element_blank(),
          plot.title = element_markdown(),
          legend.position="none")
ggsave(paste0("figures_public/", genus, "_over_time.pdf"), abund.over.time, 
       height = 6.5, width = 4)


# ### Create combined modelling dataset

# #### Combine mosquito and cattle locations

# Merge the vector and livestock data sets
mosp$popDens500m.per.bd.sds <- mos$popDens500m.per.bd.sds <- NULL
mos$Sampling <- "Vector"
mosp$Sampling <- "Cattle"
mos$ndvi8m.sds <- NA
mos$NRainDays8m.log10.sds <- NA
setdiff(names(mosp), names(mos))
setdiff(names(mos), names(mosp))
all(sort(names(mosp)) == sort(names(mos)))
mos.all <- rbind(mos, mosp[, names(mos)])
mos.all$Sampling <- factor(mos.all$Sampling)
table(mos.all$Sampling)

# Make quadratic NDVI terms
mos.all$ndvi.sds2 <- mos.all$ndvi.sds^2
mos.all$ndvi8m.sds2 <- mos.all$ndvi8m.sds^2
covariates <- c(covariates, `NDVI²` = "ndvi.sds2")

# Make quadratic rain terms
mos.all$NRainDays.log10.sds2 <- mos.all$NRainDays.log10.sds^2
covariates <- c(covariates, `n rainy days²` = "NRainDays.log10.sds2")

# Create covariates representing interactions between irrigated cropland 
# and sampling round, rain & NDVI
mos.all$IrrigatedXRain <- mos.all$proportionIrrigatedCrop.sds * mos.all$NRainDays.log10.sds

mos.all$IrrigatedXNDVI <- mos.all$proportionIrrigatedCrop.sds * mos.all$ndvi.sds
mos.all$IrrigatedXNDVI2 <- mos.all$proportionIrrigatedCrop.sds * mos.all$ndvi.sds2

mos.all$IrrigatedXAug23 <- mos.all$proportionIrrigatedCrop.sds * (mos.all$samplingRound == "Aug23")
mos.all$IrrigatedXMay24 <- mos.all$proportionIrrigatedCrop.sds * (mos.all$samplingRound == "May24")

# Add the interaction variables to the covariates list
covariates <- 
  c(covariates, 
    `irrigated:rain` = "IrrigatedXRain",
    `irrigated:NDVI` = "IrrigatedXNDVI",
    `irrigated:NDVI²` = "IrrigatedXNDVI2",
    `irrigated:Aug23` = "IrrigatedXAug23",
    `irrigated:May24` = "IrrigatedXMay24")
head(mos.all[, covariates])

# Check for extreme correlations
hist(cor(st_drop_geometry(mos.all[, covariates]), use = "pairwise.complete"),
     breaks = 1000)

# Show the range of points in km
apply(apply(st_coordinates(mos.all), 2, range), 2, diff)

# ### Exploratory spatial visualisation

# #### Sampling locations map

# Make a map of the trap locations. 
trap.map <- 
  ggplot() +
  annotation_map_tile(type = "cartolight", # available types: rosm::osm.types()
                      zoomin = 0) +  # You can change the tile type
  geom_sf(data = mos.all, 
          aes(colour = Sampling, shape = TimesSampled),
          inherit.aes = FALSE,
          alpha = 0.6,
          size = 2) +
  annotation_scale(location = "bl", width_hint = 0.25, height = unit(0.1, "cm")) +
  theme_minimal()
trap.map

# #### Observed abundance maps

# Make a map showing mean vector abundance per village at each sampling round
bb <- st_bbox(mos, crs = crs.km) * 1000
dx <- 0.05 * (bb["xmax"] - bb["xmin"])
dy <- 0.05 * (bb["ymax"] - bb["ymin"])
bb_fixed <- c(
  bb["xmin"] - dx,
  bb["xmax"] + dx,
  bb["ymin"] - dy,
  bb["ymax"] + dy)
vector.map.palette <- rev(sequential_hcl(100, "Purple-Yellow"))
vector.map.breaks <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
vector.map.list <- 
  lapply(levels(mos$samplingRound), function(s) {
    moss <- mos[mos$samplingRound == s, ]
    mossv <- 
      summarise(group_by(moss, Village), 
                geometry = st_centroid(st_combine(geometry)), 
                Total.Females = mean(Total.Females),
                Class = unique(Class),
                .groups = "drop")
    stopifnot(max(vector.map.breaks) > max(mossv$Total.Females))
    vector.map <- 
      ggplot() +
      annotation_map_tile(type = "cartolight", # available types: rosm::osm.types()
                          zoomin = 0) +  
      geom_sf(data = mossv,
              aes(shape = Class),
              colour = "black",
              size = 2.3,
              show.legend = FALSE) +
      geom_sf(data = mossv,
              aes(colour = Total.Females,
                  shape = Class),
              size = 2,
              alpha = 1) +
      coord_sf(xlim = c(bb_fixed["xmin"], bb_fixed["xmax"]),
               ylim = c(bb_fixed["ymin"], bb_fixed["ymax"]),
               expand = FALSE) +
      scale_colour_gradientn(
        colours = vector.map.palette,
        trans = "log1p",
        breaks = vector.map.breaks,
        limits = range(vector.map.breaks),
        name = "Mean count") +
      annotation_scale(location = "tr", width_hint = 0.25, height = unit(0.1, "cm")) +
      labs(tag = LETTERS[match(s, levels(mos$samplingRound)) * 2 - (genus == "Culex")]) +
      guides(
        size   = guide_legend(order = 1),
        colour = guide_colourbar(order = 2)) +
      theme_minimal()
    
    # Make a histogram showing the distribution of village mean vector count
    hist.inset <-
      ggplot(mossv,
             aes(x = Total.Females)) +
      geom_histogram(
        aes(fill = after_stat(x)),
        #bins = 15,
        breaks = vector.map.breaks,
        colour = "black",
        linewidth = 0.2) +
      scale_fill_gradientn(
        colours = vector.map.palette,
        trans = "log1p",
        guide = "none") +
      scale_x_continuous(
        trans = "log1p",
        breaks = vector.map.breaks,
        limits = range(vector.map.breaks)) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())
    
    # Add inset histogram
    vector.map.with.hist <-
      vector.map +
      inset_element(
        hist.inset,
        left   = 0.3,
        bottom = 0.15,
        right  = 0.7,
        top    = 0.55)
    vector.map.with.hist
  })
names(vector.map.list) <- levels(mos$samplingRound)

wrap_plots(vector.map.list, ncol = 1) + 
  plot_layout(guides = 'collect')

# Save the plot R objects so that they can be loaded 
# in the same script and combined
saveRDS(vector.map.list, paste0("figures_public/", genus, "_obs_count_mean.rds"))


# ## MODEL-FITTING ----

# ### Priors on fixed effects:  
# 
# See ?control.fixed. Default prior on the intercept is 
# N(mu=0, tau=0.001) (a change from INLA, which had an improper prior, to inlabru).
# This is also the default prior for the non-intercept fixed effects: N(mu=0, tau=0.001).
# Use these priors. For continuous covariates, the prior could be informative 
# depending on the scale, therefore it made sense to scale continuous covariates to SD=1.
# Centering is also advisable as this stops the intercept having extreme values
# that might stray outside the plausible prior range.
# 
# See also inla.set.control.fixed.default(). Available priors: names(inla.models()$prior).  
#   

# ### Priors for IID random effects.   
# 
# Quoting from
# https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixed.html#independent-random-effects-iid-model
# "The internal parameterization of the hyperparameter of this model
#  is \theta = log(\tau), which by default is assigned a log-Gamma distribution 
#  with parameters 1 and 0.00005. In other words, the prior of \tau is a 
#  Gamma distribution." 
#  Alternative prior choices can be set using f(..., hyper = list(...))
# -- see docs for details. https://www.r-inla.org/documentation  
# See also https://tutorials.inbo.be/tutorials/r_inla/random_intercept.pdf
#   

# The default prior appears to push the variance strongly towards 0. Try less restrictive 
# alternatives, but be careful not to penalise variances close to zero:
prior.shape <- 0.5 # default is 1
prior.rate <- 0.01 # default is 0.00005

# Visualise the prior distribution
prior.sample.prec <- rgamma(10000, prior.shape, prior.rate)
prior.sample.prec[prior.sample.prec < 0.001] <- 0.001 # truncate the sample for ease of plotting
prior.sample.var <- 1/prior.sample.prec
prior.sample.sd <- sqrt(prior.sample.var)

# Calculate quantiles to show that a broad range of variances is allowed
prop.table(table(cut(prior.sample.var, c(0, 0.01, 0.05, 4, 10, Inf))))
quantile(prior.sample.var, c(0.025, 0.5, 0.975))

# Plot the IID random effect priors
par(mfrow = c(3, 1))
hist(prior.sample.prec, xlab = "Precision", breaks = round(max(prior.sample.prec), -1),
     freq = FALSE)
lines(0:ceiling(max(prior.sample.prec)), 
      dgamma(0:ceiling(max(prior.sample.prec)), prior.shape, prior.rate), col = "pink")
abline(v = median(prior.sample.prec), col = "red") 
hist(prior.sample.var, xlab = "Variance", xlim = c(0, 50), 
     breaks = round(max(prior.sample.var), -1))
abline(v = median(prior.sample.var), col = "red") 
hist(prior.sample.sd, xlab = "SD", xlim = c(0, 3), breaks = length(prior.sample.sd))
abline(v = median(prior.sample.sd), col = "red") 
rm(prior.sample.prec, prior.sample.var, prior.sample.sd)
par(mfrow = c(1, 1))

# These priors appear realistic, so use them
iid.prec.prior <-
  list(prec = list(prior = "loggamma", param = c(prior.shape, prior.rate)))

# ### Prior for the negative binomial dispersion parameter
#    
# gamma(1, 0.2) gives reasonable flexibility to 
# have quite overdispersed counts (low theta) or nearer to Poisson counts (high theta):
prior.nb <- list(hyper = list(theta = list(prior = "loggamma", param = c(1, 0.2))))


# ### Setting up the mesh

# max.edge[1] defines the largest allowed triangle edge length in the main area  
# max.edge[2] defines the largest allowed triangle edge length in the edge area  
set.seed(801143597) # so that we always get the same mesh
hull <- fm_extensions(mos.all[mos.all$Sampling == "Vector", ], convex = c(25, 50))
Mesh <- fm_mesh_2d(mos.all[mos.all$Sampling == "Vector", ], max.edge = c(10, 40), cutoff = 0.01, boundary = hull)

# Plot the mesh
mesh.sfc <- fm_as_sfc(Mesh)
st_crs(mesh.sfc) <- st_crs(mos.all[mos.all$Sampling == "Vector", ])
ggplot(data = mos.all[mos.all$Sampling == "Vector", ]) +
  geom_sf(data = mesh.sfc) +
  geom_sf(aes(color = Sampling, shape = Sampling), alpha = 0.1, size = 3) +
  theme_minimal() 

# inla.spde2.matern() makes a Matern SPDE model object  
# The SPDE allows us to approximate in a computationally efficient way 
# a continuous Matern correlation structure with a mesh.
spde <- inla.spde2.pcmatern(mesh = Mesh, 
                            prior.range = c(1, 0.025), 
                            prior.sigma = c(1, 0.1),
                            constr = TRUE)

# Priors on the Matern correlation:  
# prior.range = c(range0, Prange)
# implies that P(range < range0) = Prange
# e.g. prior.range = c(1, 0.025) implies that we think there's 
# a 2.5% chance that the spatial range is < 1 km. 
# The spatial range is defined as the distance at which correlation
# between two locations declines to 0.1. This is quite a subjective and probably 
# influential decision, so we should do a sensitivity analysis.
# prior.sigma = c(1, 0.1)
# implies that P(sigma > sigma0) = Psigma
# e.g. prior.sigma = c(1, 0.1) implies that we think there's 
# a 10% chance that the SD of the spatial field is > 1.  
#   


# ### Making the A matrix 

# Fit a spatial model with a time-constant spatial field (i.e. the spatial correlation
# structure doesn't change over the three sampling points).  

# Name the A matrix A.tc for time-constant
A.tc <- inla.spde.make.A(mesh = Mesh, loc = mos.all[mos.all$Sampling == "Vector", ])
dim(A.tc)

# Creating the spatial field
field.tc <- inla.spde.make.index(
  name    = 'field', 
  n.spde  = spde$n.spde)  

# Allow the spatial field to vary by sampling round. This gave a better
# fitting model (lower WAIC).   

# Name the A matrix A.tv for time-varying
Groups <- "samplingRound"
NGroups <- length(unique(mos[[Groups]])) 
A.tv <- inla.spde.make.A(mesh = Mesh,
                         loc = mos.all[mos.all$Sampling == "Vector", ],
                         group = as.numeric(mos.all[mos.all$Sampling == "Vector", ][[Groups]]),
                         n.group = NGroups)

dim(A.tv)
ncol(A.tv) / NGroups

# We are now creating NGroups spatial fields
field.tv <- inla.spde.make.index(
  name    = 'field', 
  n.spde  = spde$n.spde,
  n.group = NGroups)  

# ### Making the model matrices for the fixed effects  

# Make a null model and model matrices to answer our research questions.  
Models <- list()
Models$Null <- ""
Models$E1 <- 
  "elevation500m.sds + ndvi.sds + NRainDays.log10.sds"
Models$E2 <- 
  "elevation500m.sds + ndvi.sds + ndvi.sds2 + NRainDays.log10.sds"
Models$E3 <- 
  "elevation500m.sds + ndvi.sds + NRainDays.log10.sds + NRainDays.log10.sds2"
Models$E4 <- 
  "elevation500m.sds + ndvi.sds + ndvi.sds2 + NRainDays.log10.sds + NRainDays.log10.sds2"
Models$E2B <-
  paste(Models$E2, 
        "+ buildingDensity500m.log10.sds")
Models$E2P <-
  paste(Models$E2, 
        "+ popDens500m.log10.sds")
Models$E2I <-
  paste(Models$E2, 
        "+ proportionIrrigatedCrop.sds + IrrigatedXRain")
Models$E2BP <-
  paste(Models$E2, 
        "+ buildingDensity500m.log10.sds + popDens500m.log10.sds")
Models$E2BI <-
  paste(Models$E2, 
        "+ buildingDensity500m.log10.sds + proportionIrrigatedCrop.sds + IrrigatedXRain")
Models$E2PI <-
  paste(Models$E2, 
        "+ popDens500m.log10.sds + proportionIrrigatedCrop.sds + IrrigatedXRain")
Models$E2BPI <-
  paste(Models$E2, 
        "+ buildingDensity500m.log10.sds + popDens500m.log10.sds + proportionIrrigatedCrop.sds + IrrigatedXRain")

# Create a model matrix for each fixed effects model
X.list <- 
  lapply(Models, function(m) {
    if(m == "") return(NA) else
      as.data.frame(model.matrix(as.formula(paste0(" ~ 1 + ", m)), 
                                 data = mos.all[mos.all$Sampling == "Vector", ]))[, -1]
  })
sapply(X.list, dim)
sapply(X.list, function(x) sum(is.na(x)))
unique(unlist(lapply(X.list, colnames)))

# ### Fit a negative binomial GLMM for each model matrix (takes around 10 min)
# 
# Each model matrix represents a different fixed effects model.
# Also fit a zero-inflated model only for the null fixed effects model.
# (Takes 5-10 min)
inla.opt <- 
  list(
    num.threads = "1:1",
    inla.mode = "compact",
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                           config = TRUE,
                           return.marginals.predictor = TRUE,
                           openmp.strategy = NULL),
    control.family = prior.nb)


# Set random seed for reproducibility
set.seed(420781748, kind = "L'Ecuyer-CMRG") # robust parallel-safe RNG

# Fit each model
Model.fits <-
  mclapply(names(Models), function(modname) {
    
    start.model <- Sys.time()
    message_parallel(paste("Starting model", modname, "at", format(start.model, "%H:%M")))
    
    # Specify the model components
    # The sampling round random effect has a fairly narrow prior 
    # to prevent extreme variances given that there are only three rounds.
    # list(prior = "pc.prec", param = c(5, 0.01)) means that there is a 
    # 1% chance that the sampling round SD exceeds 5. 
    components <- ~ Intercept(1) +
      field.tc(geometry, model = spde) +
      field.tv(geometry, model = spde, group = as.numeric(samplingRound),
               control.group = list(model = "iid")) +
      field.tve(geometry, model = spde, group = as.numeric(samplingRound),
                control.group = list(model = "exchangeable")) +
      Village(Village, model = "iid", constr = TRUE, hyper = iid.prec.prior) +
      VillageTV(VillageTV, model = "iid", constr = TRUE, hyper = iid.prec.prior) +
      samplingRound(samplingRound, model = "iid", constr = TRUE,
                    hyper = list(prec = list(prior = "pc.prec", param = c(5, 0.01)))) + 
      X(main = X.list[[modname]], model = "fixed")
    
    # Specify the observation model
    obs <- 
      if(modname == "Null") {
        bru_obs(
          formula = Total.Females ~ Intercept + field.tve + samplingRound,
          family = "nbinomial",
          data = mos.all[mos.all$Sampling == "Vector", ])
      } else {
        bru_obs(
          formula = Total.Females ~ Intercept + field.tve + samplingRound + X,
          family = "nbinomial",
          data = mos.all[mos.all$Sampling == "Vector", ])
      }
    
    # Fit the models
    fit <- 
      bru_rerun(bru(
        components = components,
        obs = obs,
        options = inla.opt))
    fit$waic$waic
    
    # Fit alternative random effects specification to the Null fixed effects model
    if(modname %in% c("Null")) {
      # Only IID random effects
      obs.nonspatial.iid <- 
        bru_obs(
          formula = Total.Females ~ Intercept + Village + VillageTV + samplingRound,
          family = "nbinomial",
          data = mos.all[mos.all$Sampling == "Vector", ])
      RE1 <- 
        bru_rerun(bru(
          components = components,
          obs = obs.nonspatial.iid,
          options = inla.opt))
      
      # Time-constant spatial random effect 
      obs.spatial.tc <- 
        bru_obs(
          formula = Total.Females ~ Intercept + field.tc + samplingRound,
          family = "nbinomial",
          data = mos.all[mos.all$Sampling == "Vector", ])
      RE2 <- 
        bru_rerun(bru(
          components = components,
          obs = obs.spatial.tc,
          options = inla.opt))
      
      # Time-constant spatial random effect
      # and with IID random effects
      obs.spatial.iid.tc <- 
        bru_obs(
          formula = Total.Females ~ Intercept + Village + VillageTV + field.tc + samplingRound,
          family = "nbinomial",
          data = mos.all[mos.all$Sampling == "Vector", ])
      RE3 <- 
        bru_rerun(bru(
          components = components,
          obs = obs.spatial.iid.tc,
          options = inla.opt))
      
      # Time-varying spatial random effect with correlation between the three fields
      obs.spatial.iid.tve <- 
        bru_obs(
          formula = Total.Females ~ Intercept + Village + VillageTV + field.tve + samplingRound,
          family = "nbinomial",
          data = mos.all[mos.all$Sampling == "Vector", ])
      RE5 <- 
        bru_rerun(bru(
          components = components,
          obs = obs.spatial.iid.tve,
          options = inla.opt))
      
      # Time-varying spatial random effect with correlated fields
      # (already fitted above)
      RE4 <- fit
      
    } else RE1 <- RE2 <- RE3 <- RE4 <- RE5 <- NULL
    
    attr(fit, "re.models") <- list(
      RE1 = RE1,
      RE2 = RE2,
      RE3 = RE3,
      RE4 = RE4,
      RE5 = RE5)
    
    # Fit a ZI NB alternative to the Null fixed effects NB model
    if(modname == "Null") {
      obs.zi <- 
        bru_obs(
          formula = Total.Females ~ Intercept + field.tve + samplingRound,
          family = "zeroinflatednbinomial1",
          data = mos.all[mos.all$Sampling == "Vector", ])
      fit.zi <- 
        bru_rerun(bru(
          components = components,
          obs = obs.zi,
          options = inla.opt))
    } else fit.zi <- NULL
    
    attr(fit, "fit.zi") <- fit.zi
    
    message_parallel(paste("Finished model", modname, "after", Sys.time() - start.model))
    
    return(fit)
  }, mc.cores = detectCores() - 1, mc.set.seed = FALSE)
names(Model.fits) <- names(Models)

# Convergence plot (comment out for speed)
#lapply(Model.fits, bru_convergence_plot)

# Check all models converged
all(sapply(Model.fits, "[[", "ok"))

# Check if any models have suspiciously large SDs for fixed effects
sapply(Model.fits[-1], function(m) max(m$summary.random$X$sd, na.rm = TRUE))

# Check for unusually large KLD in the random and fixed effects
sapply(Model.fits, function(m) max(unlist(lapply(m$summary.random, "[[", "kld"))))

# Look for outliers in marginal log-likelihood
range(sapply(Model.fits, function(m) m$mlik[1, 1]))

# WAIC
IC.inla(Model.fits, crit = "waic")

# ## MODEL VALIDATION ----

# ### Assessing the validity of the model assumptions

# ### Use residuals of the null model to assess goodness-of-fit

# Calculate the fixed effects predictions and the spatial random effect
fit0.linpred <- Model.fits$Null$summary.fixed$mean

# Calculate the spatial random effect predictions
spatial.re.mean <- matrix(Model.fits$Null$summary.random$field.tve$mean, nrow = spde$n.spde, ncol = NGroups)
spatial.predictions <- 
  apply(spatial.re.mean, 2, function(field) fm_evaluate(fm_evaluator(Mesh, loc = mos.all[mos.all$Sampling == "Vector", ]), field))

# Combine the predictions from the  fixed and random effects
fit0.pred.cond <- 
  fit0.linpred +
  spatial.predictions[cbind(1:nrow(spatial.predictions), as.numeric(mos.all$samplingRound[mos.all$Sampling == "Vector"]))] #+

# Calculate expected residuals by simulation
fit0.resid.sim <- 
  sapply(1:1000, function(i) {
    fit0.sim <- 
      rnbinom(sum(!is.na(mos.all$Total.Females)), mu = exp(fit0.pred.cond), 
              size = Model.fits$Null$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", "mean"])
    sort(fit0.sim - fit0.pred.cond)
  })

# Calculate observed residuals
fit0.resid <- na.omit(mos.all[["Total.Females"]]) - fit0.pred.cond
sd(fit0.resid, na.rm = TRUE)

# Use the (simulated) expected residuals to standardise the observed 
# residuals to quantiles
unif.fit0.resid.obs <-
  sapply(sort(fit0.resid), function(x) {
    mean(x > fit0.resid.sim)
  })
unif.fit0.resid.exp <- ((1:length(unif.fit0.resid.obs)) - 0.5)/length(unif.fit0.resid.obs)

# Plot the standardised residuals. If the negative binomial distribution fits well
# they should be uniform. Assess this using a histogram, and a QQ plot of 
# observed vs expected residuals.
resid.plot.pch <- 16
resid.plot.cex <- 0.5
resid.plot.col <- alpha("blue", 0.2)
#hist(unif.fit0.resid.obs, breaks = 10, main = "Histogram of standardised residuals")
plot(unif.fit0.resid.exp, unif.fit0.resid.obs, main = "QQ plot of of standardised residuals", 
     cex = resid.plot.cex, col = resid.plot.col, pch = resid.plot.pch)
abline(0, 1, col = "red")
r.fit0 <- cor(unif.fit0.resid.exp, unif.fit0.resid.obs)
legend("topleft", legend = paste0("r=", round(r.fit0, 4)), bty = "n")
plot(jitter(rep(fit0.linpred, length(unif.fit0.resid.obs))), unif.fit0.resid.obs, 
     main = "Plot of standardised residuals\nagainst (jittered) fitted values")

# ### Posterior predictive checks

# #### Simulate replicated datasets

# Use simulated counts from the null model to assess goodness-of-fit.  
# Here we are not just assessing whether the residuals fit a negative binomial
# distribution, as above, but are simulating from the whole model, and therefore 
# including uncertainty from the whole model. The idea is to check for 
# discrepancies between the observed and simulated data that would indicate 
# that the model is a poor representation of the observed data.  
#   


#  Draw posterior samples
n.samp <- 15
samples <- generate(Model.fits$Null, newdata = mos.all[mos.all$Sampling == "Vector", ], 
                    n.samples = n.samp, seed = 53300846)

# Combine the n.samp posterior samples to produce simulated counts
simulated.y <- 
  sapply(1:length(samples), function(i) {
    
    # Spatial random field
    field.matrix <- fm_evaluate(fm_evaluator(Mesh, loc = mos.all[mos.all$Sampling == "Vector", ]), 
                                field = matrix(samples[[i]]$field.tve, ncol = NGroups))
    field.values <- field.matrix[cbind(1:nrow(field.matrix), as.numeric(mos.all$samplingRound[mos.all$Sampling == "Vector"]))]
    
    # Add spatial RE to fixed effects and other random effects
    eta <- samples[[i]]$Intercept + field.values# +
    
    # Finally simulate counts from a negative binomial distribution
    y <- rnbinom(length(eta), mu = exp(eta), 
                 size = samples[[i]]$`size_for_the_nbinomial_observations_1/overdispersion_`)
    return(y)
  })

# Tabulate the simulated counts beside the observed counts
simtab <- cbind(Observed = mos.all$Total.Females[mos.all$Sampling == "Vector"], 
                Simulated = simulated.y)
dim(simtab) # should be n_obs x (1 + n.samp)
# Bin the very rare higher counts to facilitate visualisation
fac.lev <- 
  c(0:20, 
    seq(21, 101, 10),
    seq(201, 1001, 100), 
    seq(2001, max(2001, max(mos.all$Total.Females, simtab, na.rm = TRUE)), 1000), 
    Inf) - 0.5
names(fac.lev) <- paste(fac.lev[-length(fac.lev)] + 0.5, fac.lev[-1] - 0.5, sep = "-")
names(fac.lev)[1:21] <- 0:20
simtab.plotdata <-
  t(apply(simtab, 2, function(x) {
    table(cut(x, fac.lev, labels = names(fac.lev)[-length(fac.lev)]))
  }))
dim(simtab.plotdata)


# Barplots comparing distributions of observed counts and posterior samples
old.par <- par(mfrow = c(5, 1), mar = c(5.1, 2, 0, 2))
sapply(1:ceiling(ncol(simtab.plotdata)/11), function(i) {
  idx <- (i*11 - 10):min(ncol(simtab.plotdata), (i*11))
  if(sum(simtab.plotdata[, idx]) > 0) {
    barplot(simtab.plotdata[, idx],
            beside = TRUE,
            cex.names = 0.8,
            las = 2,
            col = c("grey", "red")[1 + (colnames(simtab) == "Observed")])
  }
  NULL
})
par(old.par)

# Boxplots to compare observed and simulated counts broken down by sampling round
old.par <- par(mfrow = c(4, 4), mar = c(2, 2, 2, 2))
plot(log10(1 + Total.Females) ~ samplingRound, 
     data = mos.all[mos.all$Sampling == "Vector", ],
     main = "Observed", ylim = c(0, log10(1 + max(Total.Females, simulated.y))),
     col = "red")
lapply(1:ncol(simulated.y), function(i) {
  plot(log10(1 + simulated.y[, i]) ~ samplingRound, 
       data = mos.all[mos.all$Sampling == "Vector", ],
       main = "Simulated", ylim = c(0, log10(1 + max(Total.Females, simulated.y))))
})
par(old.par)


# Compare the posterior samples to the observed counts using summaries
# (posterior predictive checks)
old.par <- par(mfrow = c(2, 2))
hist(apply(simulated.y, 2, sd), main = "Histogram of SD of posterior samples")
abline(v = sd(mos.all$Total.Females, na.rm = TRUE), col = 2)
legend("topright", lty = 1, col = 2, legend = "Observed SD", bty = "n")
hist(apply(simulated.y, 2, mean), main = "Histogram of mean of posterior samples")
abline(v = mean(mos.all$Total.Females, na.rm = TRUE), col = 2)
legend("topright", lty = 1, col = 2, legend = "Observed mean", bty = "n")
hist(apply(simulated.y, 2, function(x) diff(range(x))), 
     main = "Histogram of range of posterior samples")
abline(v = diff(range(mos.all$Total.Females, na.rm = TRUE)), col = 2)
legend("topright", lty = 1, col = 2, legend = "Observed range", bty = "n")
hist(apply(simulated.y, 2, function(x) acf(x, plot = FALSE)$acf[2]), 
     main = "Histogram of ACF(1) of posterior samples")
abline(v = acf(na.omit(mos.all$Total.Females), plot = FALSE)$acf[2], col = 2)
legend("topright", lty = 1, col = 2, legend = "Observed ACF(1)", bty = "n")
par(old.par)


# Plot posterior samples against observed counts
plot(1 + mos.all$Total.Females, 1 + mos.all$Total.Females, 
     type = "n",
     xlab = "Observed count + 1",
     ylab = "Simulated count + 1",
     xlim = 1 + range(mos.all$Total.Females, na.rm = TRUE),
     ylim = 1 + range(simulated.y),
     log = "xy")
apply(simulated.y, 2, function(y) points(1 + na.omit(mos.all$Total.Females), 1 + y))
abline(0, 1, col = 2)


# ### Compare the random effects specification
re.waic <- sapply(attributes(Model.fits$Null)$re.models, function(x) x$waic$waic)
formattable(data.frame(WAIC = round(re.waic, 2), '\u0394WAIC' = round(re.waic - min(re.waic), 2)))
lapply(attributes(Model.fits$Null)$re.models, function(x) x$summary.hyperpar)

# Compare (approximate) prior and posterior of rho (~1 min)
t(rbind(
  Prior = inla.na.fit(attributes(Model.fits$Null)$re.models$RE4)$summary.hyperpar["GroupRho for field.tve", ],
  Posterior = attributes(Model.fits$Null)$re.models$RE4$summary.hyperpar["GroupRho for field.tve", ]))

# ### Finally, is there evidence of zero-inflation?

# Check the ZI parameter estimate and compare ZI NB and non-ZI NB using WAIC
attributes(Model.fits$Null)$fit.zi$summary.hyperpar["zero-probability parameter for zero-inflated nbinomial_1", ]
IC.inla(list(NB = Model.fits$Null, ZINB = attributes(Model.fits$Null)$fit.zi), crit = "waic")
# Culex: The NB model has lower WAIC than ZINB, higher marginal likelihood, and the estimate 
# of ZI probability is very close to zero, so we reject ZI.  
# Aedes: The NB model has slightly (0.81, which is negligible) higher WAIC than ZINB, 
# but also higher marginal likelihood, and the estimate of ZI probability is fairly 
# close to zero (7%), so we reject ZI.

# Conclusion: The negative binomial model fits well enough for both genera.  
# 

# ### Blocked k-fold cross-validation

# #### Define spatial folds
# 
# We are using blocked cross-validation due to the non-independence problem 
# with randomly selecting held-out folds:
# Roberts et al. (2017). Cross-validation strategies for data with temporal, 
# spatial, hierarchical, or phylogenetic structure. Ecography.
# Valavi et al. (2019). blockCV: An R package for generating spatially or 
# environmentally separated folds for cross-validation of species distribution 
# models. Methods in Ecology and Evolution.   

# The aim is to choose clusters that are spatially and temporally independent.
# We do this by first estimating the spatial scale of correlation
Model.fits$Null$summary.hyperpar["Range for field.tve", "mean"]
# 12 km for Culex, 26 km for Aedes.
# A robust way to do this is by k-means clustering. 
# We want the diameter of the held-out clusters of traps to be at least double that scale.
# More folds (higher k) will give smaller clusters. 5 and 10 are common choices for k.
# I tried both, and only k=5 gave a minimum cluster diameter that was > 4 * 12 km. The code 
# for calculating cluster diameter (2 * distance from fold point centroid to nearest 
# outside point) is in the CV loop below.  

# Assign and plot clusters
cluster.k <- 5
cv <- cv_cluster(x = mos, k = cluster.k)
cv_plot(cv, mos)

# Store cluster ID in mos.all
mos.all$fold <- NA
mos.all$fold[1:nrow(cv$biomod_table)] <- 
  ((1 - cv$biomod_table) %*% 1:cluster.k)[, 1]
table(mos.all$fold)

# Make a map of the trap locations showing folds
ggplot() +
  annotation_map_tile(type = "cartolight", # available types: rosm::osm.types()
                      zoomin = 0) +  # You can change the tile type
  geom_sf(data = mos.all[!is.na(mos.all$fold), ], 
          aes(colour = factor(fold)),
          inherit.aes = FALSE,
          alpha = 0.6,
          size = 2) +
  annotation_scale(location = "bl", width_hint = 0.25, height = unit(0.1, "cm")) +
  theme_minimal()

# Set inla options for CV
inla.opt.cv <- inla.opt
inla.opt.cv$control.compute$cpo <- inla.opt.cv$control.compute$dic <-
  inla.opt.cv$control.compute$waic <- FALSE
inla.opt.cv$control.compute$return.marginals.predictor <- FALSE
inla.opt.cv$num.threads <- detectCores()
inla.opt.cv$inla.mode <- "compact"
inla.opt.cv$control.inla = list(strategy = "simplified.laplace", int.strategy = "eb", h = 0.02)

# Loop over each fold, fitting model to training folds and predicting 
# vector abundance on held-out fold
cv.results <- 
  sapply(names(Model.fits), function(fit.name) {
    
    mos.cv <- mos.all[mos.all$Sampling == "Vector", ]
    mos.cv$pred.fold <- NA
    mos.cv$fold.diameter <- NA
    
    print(fit.name)
    fit <- Model.fits[[fit.name]]
    
    # Original components
    components.cv <- as_bru_comp_list(fit)  
    
    # Loop over all k clusters
    for(i in 1:cluster.k) {
      print(i)
      print(table(!is.na(mos.cv$pred.fold)))
      
      # Calculate distance from centroid of fold to nearest point outside fold
      fold.centroid <- colMeans(st_coordinates(mos.cv[mos.cv$fold %in% i, ]))
      mos.cv$fold.diameter[mos.cv$fold %in% i] <- 
        2 * min(st_distance(st_sfc(st_point(fold.centroid), crs = crs.km), 
                            mos.cv[mos.cv$fold %in% (1:ncol(cv$biomod_table))[-i], ]))
      print(unique(mos.cv$fold.diameter[mos.cv$fold %in% i]))
      
      # Refit model to 4 of the 5 folds, leaving 1 fold out (fold i)
      inla.opt.cv$verbose <- inla.opt$verbose <- FALSE
      mos.cv.na <- mos.cv
      mos.cv.na$Total.Females[mos.cv.na$fold %in% i] <- NA
      
      # Specify the observation model
      obs.cv <- 
        if(fit.name == "Null") {
          bru_obs(
            formula = Total.Females ~ Intercept + samplingRound + field.tve,
            family = "nbinomial",
            data = mos.cv.na)
        } else {
          bru_obs(
            formula = Total.Females ~ Intercept + samplingRound + field.tve + X,
            family = "nbinomial",
            data = mos.cv.na)
        }
      
      # Fit the model
      cv.fit.start <- Sys.time()
      fit.cv <- 
        bru(
          components = components.cv,
          obs = obs.cv,
          options = inla.opt)
      cv.fit.finish <- Sys.time()
      print(cv.fit.finish - cv.fit.start)
      
      # Predict left-out fold from model fitted to other 4 folds
      pred.mean.cv <- fit.cv$summary.fixed$mean
      if(fit.name == "Null") {
        pred.mean.cv <- fit.cv$summary.fixed$mean   # intercept
      } else {
        pred.mean.cv <- 
          fit.cv$summary.fixed$mean +
          (as.matrix(X.list[[fit.name]]) %*% fit.cv$summary.random$X$mean)[which(!cv$biomod_table[, i]), 1]  
      }
      
      mos.cv$pred.fold[which(!cv$biomod_table[, i])] <- pred.mean.cv
      print(table(!is.na(mos.cv$pred.fold)))
    }
    
    tapply(mos.cv$pred.fold, mos.cv$fold, mean, na.rm = TRUE)
    
    plot(mos.cv$pred.fold, log(mos.cv$Total.Females + 0.5), main = fit.name)
    rho <- cor(mos.cv$pred.fold, mos.cv$Total.Females, use = "pairwise.complete",
               method = "spearman")
    rho.village <- 
      cor(
        tapply(mos.cv$Total.Females, mos.cv$VillageTV, mean, na.rm = TRUE),
        tapply(mos.cv$pred.fold, mos.cv$VillageTV, mean, na.rm = TRUE),
        use = "pairwise", method = "spearman")
    plot(
      tapply(log(0.5 + mos.cv$Total.Females), mos.cv$VillageTV, mean, na.rm = TRUE),
      tapply(mos.cv$pred.fold, mos.cv$VillageTV, mean, na.rm = TRUE),
      pch = unlist(tapply(mos.cv$fold, mos.cv$VillageTV, unique)), 
      main = paste0(fit.name, "\nrho=", round(rho, 2), "\nrho.village=", round(rho.village, 2)))
    
    c(rho = rho, rho.village = rho.village)
  })


# ## MODEL COMPARISON AND INFERENCE ----

# ### Tabulate and plot results for all models

# Get WAIC for all models
WAIC.tab <- cbind(Fixed_effects = unlist(Models), IC.inla(Model.fits, crit = "waic"))
for(i in 1:nrow(WAIC.tab)) {
  for(j in 1:length(covariates)) {
    WAIC.tab[i, "Fixed_effects"] <- 
      gsub(covariates[j], names(covariates)[j], WAIC.tab[i, "Fixed_effects"])
  }
  WAIC.tab[i, "Fixed_effects"] <- gsub("\\.sds", "", WAIC.tab[i, "Fixed_effects"])
}
names(WAIC.tab) <- gsub("delta", "\u0394", names(WAIC.tab))

# Add 5-fold CV results to WAIC table
WAIC.tab <- cbind(WAIC.tab, t(round(cv.results, 2))[rownames(WAIC.tab), ])

# Get variance decomposition stats for all models
Model.varcomp.list <- 
  lapply(names(Models), function(modname) {
    print(modname)
    variance.decomp(result = Model.fits[[modname]], mm = X.list[[modname]],
                    spatial.re = "field.tve",
                    fixed.est.mean = Model.fits[[modname]]$summary.random$X$mean,
                    random.effects = c("samplingRound"),
                    nb.ol.var.method = "exclude")
  }) 

Model.varcomp.tab <- t(do.call("cbind", Model.varcomp.list))
rownames(Model.varcomp.tab) <- names(Models)

# Plot the variance components for each model
old.par <- par(mar = par()$mar + c(5, 0, 0, 0))
barplot.data <-
  barplot(Model.varcomp.tab[, 1:4],
          beside = TRUE, names.arg = rep(rownames(Model.varcomp.tab), 4),
          las = 2, ylab = "Variance", 
          main = "Variance components",
          ylim = c(0, max(Model.varcomp.tab[, 1:4]) * 1.2))
text(x = apply(barplot.data, 2, mean), 
     y = rep(max(Model.varcomp.tab[, 1:4]), 4) * 1.1, 
     labels = colnames(Model.varcomp.tab)[1:4])
par(old.par)

# Plot the WAIC for each model
old.par <- par(mar = par()$mar + c(5, 0, 0, 0))
barplot(WAIC.tab[, "\u0394WAIC"], beside = TRUE, names.arg = rownames(WAIC.tab),
        las = 2, ylab = "\u0394WAIC",
        main = "Model comparison using WAIC")
par(old.par)

# Combine WAIC and variance decomposition stats in one table
Model.sum.tab <-
  cbind(
    WAIC.tab, 
    n_par = sapply(X.list, function(x) length(na.omit(x))),
    round(Model.varcomp.tab[rownames(WAIC.tab), c("Fixed_effects_var", "R2.marginal")], 2))

# ### Model comparison table
formattable(Model.sum.tab)
saveWidget(as.htmlwidget(formattable(Model.sum.tab)), 
           file = paste0("Vector_analysis/public_outputs/", genus, "_model_comparison_tab.html"),
           selfcontained = TRUE)
unlink(paste0("Vector_analysis/public_outputs/", genus, "_model_comparison_tab_files"), 
       recursive = TRUE)

# ### Fixed-effect inference

# #### Forest plots of fixed effects estimates
forest.plots <-
  lapply(names(Model.fits), function(modname) {
    fit <- Model.fits[[modname]]
    fit.summary.fixed <- 
      rbind(cbind(ID = "Intercept", fit$summary.fixed), fit$summary.random$X)
    rm(fit)
    
    # make variable names prettier
    fit.summary.fixed$ID[-1] <- names(covariates)[match(fit.summary.fixed$ID, covariates)][-1]

    fit.summary.fixed$ID <- factor(fit.summary.fixed$ID, fit.summary.fixed$ID)
    rownames(fit.summary.fixed) <- fit.summary.fixed$ID
    
    # exponentiate estimates
    fit.summary.fixed$estimate <- exp(fit.summary.fixed$mean)
    fit.summary.fixed$lower <- exp(fit.summary.fixed$`0.025quant`)
    fit.summary.fixed$upper <- exp(fit.summary.fixed$`0.975quant`)
    
    # Export the fixed effect estimates.
    write.csv(fit.summary.fixed, 
              file = paste0("Vector_analysis/public_outputs/fixed_effects_estimates/", 
                            genus, ".", modname, ".fixed.csv"), 
              row.names = FALSE)
    
    # Plot
    out <- 
      ggplot(fit.summary.fixed, aes(y = ID, x = estimate)) +
      geom_pointrange(aes(xmin = lower, xmax = upper)) +
      geom_vline(xintercept = 1, linetype = 2) +
      scale_y_discrete(limits = rev(levels(fit.summary.fixed$ID))) +
      scale_x_continuous(breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
                         transform = "log") +
      theme_minimal() +
      labs(y = "", x = "Estimate \u00B1 95% CrI",
           title = modname)
    out
  })
names(forest.plots) <- names(Model.fits)
ggsave(paste0("figures_public/", genus, "_forest_plots.pdf"), wrap_plots(forest.plots, ncol = 2), 
       height = 11, width = 11)

# View the fixed effects 95% CrIs from each model
lapply(Model.fits[names(Model.fits) != "Null"], 
       function(m) 
         cbind(ID = m$summary.random$X$ID, 
               exp(m$summary.random$X[, c("0.025quant", "0.975quant")])))

# Make fixed effects-based predictions from all models based on the best environmental model
# to see how closely correlated they are
all.mod.pred <- 
  sapply(grep("E2", rownames(Model.sum.tab), value = TRUE)[-1], function(modname) { 
    as.matrix(X.list[[modname]]) %*% Model.fits[[modname]]$summary.random$X$mean
  })

dim(all.mod.pred)
plot(as.data.frame(all.mod.pred))
min(cor(all.mod.pred))
round(cor(all.mod.pred), 2)

# ### Spatial inference

# #### Posterior Matern correlation functions
posterior.matern.plots <-
  lapply(names(Model.fits), function(modname) {
    fit <- Model.fits[[modname]]
    post.dat <- spde.posterior(fit, "field.tve", what = "matern.correlation")
    spat.range <- post.dat[which.min(abs(post.dat$q0.5 - 0.1)), c("x", "q0.5")]
    ggplot(data = post.dat) +
      geom_ribbon(aes(x = x, ymin = `q0.025`, ymax = `q0.975`), fill = "grey") +
      geom_line(aes(x = x, y = q0.5)) +
      geom_path(data = data.frame(x = c(0, spat.range$x, spat.range$x),
                                  y = c(0.1, 0.1, 0)), 
                mapping = aes(x = x, y = y), linetype = 2) +
      xlab("Distance (km)") +
      ylab("Matern correlation") +
      ggtitle(modname) + 
      theme_minimal()
  })

# Extract the maximum x.range[2] and apply it to all plots
matern.x.lim <-
  sapply(posterior.matern.plots, 
         function(p) ggplot_build(p)$layout$panel_params[[1]]$x.range[2])
posterior.matern.plots <-
  lapply(posterior.matern.plots, function(p) p + xlim(0, max(matern.x.lim)))

# Export Matern plots
ggsave(filename = paste0("figures_public/", genus, "_matern_plots.pdf"),
       plot = wrap_plots(posterior.matern.plots), 
       height = 7, width = 11)

# Choose best model, based on WAIC and 5-fold CV
if(genus == "Culex") Best <- "E2I"
if(genus == "Aedes") Best <- "E2B"

# Visualise the relationship between NDVI and vector abundance
NDVIsds <- seq(min(mos$ndvi.sds), max(mos$ndvi.sds), 0.01)
plot(NDVIsds * attr(mos$ndvi.sds, "scaled:scale") + attr(mos$ndvi.sds, "scaled:center"), 
     NDVIsds * Model.fits[[Best]]$summary.random$X[Model.fits[[Best]]$summary.random$X$ID == "ndvi.sds", "mean"] + 
       (NDVIsds^2) * Model.fits[[Best]]$summary.random$X[Model.fits[[Best]]$summary.random$X$ID == "ndvi.sds2", "mean"], 
     type = "l", xlab = "NDVI", ylab = "log predicted vector abundance")
title(genus)

# Visualise the interaction between irrigation and rain (if present in best model)
X.new <- 
  expand.grid(elevation500m.sds = 0, ndvi.sds = 0, ndvi.sds2 = 0, 
              buildingDensity500m.log10.sds = 0, popDens500m.log10.sds = 0,
              proportionIrrigatedCrop.sds = c(min(mos$proportionIrrigatedCrop.sds),
                                              mean(range(mos$proportionIrrigatedCrop.sds)),
                                              max(mos$proportionIrrigatedCrop.sds)),
              NRainDays.log10.sds = seq(min(mos$NRainDays.log10.sds), max(mos$NRainDays.log10.sds), length.out = 50))
X.new$IrrigatedXRain <- X.new$NRainDays.log10.sds * X.new$proportionIrrigatedCrop.sds
X.new$pred.mean <- 
  Model.fits[[Best]]$summary.fixed$mean +
  (as.matrix(X.new[, Model.fits[[Best]]$summary.random$X$ID]) %*% 
     Model.fits[[Best]]$summary.random$X$mean)[, 1]
X.new$proportionIrrigatedCrop.sds.fac <- factor(X.new$proportionIrrigatedCrop.sds)
levels(X.new$proportionIrrigatedCrop.sds.fac) <- c("Low", "Medium", "High")
ggplot(data = X.new, mapping = aes(x = NRainDays.log10.sds, y = pred.mean, 
                                   group = proportionIrrigatedCrop.sds.fac,
                                   linetype = proportionIrrigatedCrop.sds.fac)) +
  geom_line()

# ## PREDICT AND EXPORT VECTOR ABUNDANCE ----

# ### Predict quarterly abundance

# ### Generate quarterly predictions

# Predict vector abundance for all sampling locations only from fixed effects.
# Problem: The relationship between the dynamic variables (NDVI and rain days)
# has been learned from the three months preceding vector trapping, over a 500m
# radius. We want to use this relationship to estimate past vector exposure
# over the lifespans of cattle ranging from 1-15 years.
# Solution: Estimate log vector abundance using the coefficients from the vector 
# model for each quarter covering the 15-year period during which the oldest
# animal in the sample lived. Then average over these quarterly abundance estimates
# for the lifespan of each animal.


# Make model matrix for cattle
X <-
  model.matrix(as.formula(paste0(" ~ 1 + ", Models[[Best]])),
               data = mos.all)
dim(X)
rownames(X) <- mos.all$newID

# ### Estimate lifetime exposure

# The prediction is linear in the quarterly terms, so predicting from the
# lifetime averages equals averaging the quarterly predictions.
# Replace placeholder values for NDVI, NDVI^2 and rain days with lifetime averages
X[match(cattle$newID, mos.all$newID), "ndvi.sds"] <- cattle$ndvi.sds
X[match(cattle$newID, mos.all$newID), "ndvi.sds2"] <- cattle$ndvi.sds2
X[match(cattle$newID, mos.all$newID), "NRainDays.log10.sds"] <- cattle$NRainDays.log10.sds

# Extract fixed effects posterior means from best model
B  <- c(Model.fits[[Best]]$summary.fixed$mean, Model.fits[[Best]]$summary.random$X$mean)

# Lifetime mean predicted log vector abundance (traps filled in below)
mos.all[, paste0("mean", genus)] <-
  ifelse(is.na(mos.all$birthyear), NA, (X %*% B)[, 1])

mos.all[[paste0("mean", genus)]][mos.all$Sampling == "Vector"] <- 
  predict(Model.fits[[Best]], newdata = mos.all[mos.all$Sampling == "Vector", ],
          formula = ~ Intercept + X, n.samples = 10000)$mean
mos.all$PredictedCount <- 
  exp(mos.all[[paste0("mean", genus)]] + sum(Model.varcomp.tab[Best, c("RE_var.Spatial", "RE_var.samplingRound")])/2)


# Predict vector abundance for cattle incidence data
# Make model matrix for cattle incidence data
X.inc <- 
  model.matrix(as.formula(paste0(" ~ 1 + ", Models[[Best]])), 
               data = cattle.inc)
dim(X.inc)
rownames(X.inc) <- cattle.inc$newID
# multiply model matrix by coefficients to produce predicted log vector abundance
cattle.inc[, paste0("mean", genus)] <-
  X.inc %*% c(Model.fits[[Best]]$summary.fixed$mean, Model.fits[[Best]]$summary.random$X$mean)


# Write cattle incidence data to CSV
write.csv(cattle.inc, 
          paste0("data_public/cattleIncidenceDataWith", genus, "Pred.csv"), 
          row.names = FALSE)

# #### Compare observed and predicted abundance

# Plot observed against predicted count
plot.default(1 + st_drop_geometry(mos.all[mos.all$Sampling == "Vector", c("Total.Females", "PredictedCount")]), 
             log = "xy")
abline(0, 1)
boxplot(1 + st_drop_geometry(mos.all[mos.all$Sampling == "Vector", c("Total.Females", "PredictedCount")]),
        log = "y")

# Note that we would expect the predicted counts to be strongly shrunk towards the mean
# because we are averaging over all the random effects, which are the spatial field and
# sampling round:
Model.fits[[Best]]$summary.hyperpar

# However, the means should be not too different: 
sapply(st_drop_geometry(mos.all[mos.all$Sampling == "Vector", c("Total.Females", "PredictedCount")]), mean)

# ### Predicted vector abundance maps

# Make maps of the vector count predictions at the cattle sampling locations
pred.map.points <- mos.all[mos.all$Sampling == "Cattle", c("Class", "PredictedCount")]
pred.map.points$geometry.fac <- as.character(pred.map.points$geometry)
pred.map.points$PredictedCount.mean <- 
  round(tapply(pred.map.points$PredictedCount, pred.map.points$geometry.fac, mean)[as.character(pred.map.points$geometry.fac)], 1)
pred.map.points$PredictedCount <- NULL
dim(unique(pred.map.points))
stopifnot(max(vector.map.breaks) > max(pred.map.points$PredictedCount.mean))
pred.map <- 
  ggplot() +
  annotation_map_tile(type = "cartolight", # available types: rosm::osm.types()
                      zoomin = 0) +  
  geom_sf(data = unique(pred.map.points),
          aes(shape = Class),
          colour = "black",
          size = 2.3,
          show.legend = FALSE) +
  geom_sf(data = unique(pred.map.points),
          aes(colour = PredictedCount.mean,
              shape = Class),
          size = 2,
          alpha = 1) +
  scale_colour_gradientn(
    colours = vector.map.palette,
    trans = "log1p",
    breaks = vector.map.breaks,
    limits = range(vector.map.breaks),
    name = "Predicted count") +
  annotation_scale(location = "br", width_hint = 0.25, height = unit(0.1, "cm")) +
  guides(
    size   = guide_legend(order = 1),
    colour = guide_colourbar(order = 2)) +
  ggtitle(paste("Predicted", genus,"abundance at cattle sampling locations")) +
  theme_minimal()

# Make a histogram showing the distribution of 
hist.inset <-
  ggplot(unique(pred.map.points),
         aes(x = PredictedCount.mean)) +
  geom_histogram(
    aes(fill = after_stat(x)),
    breaks = vector.map.breaks,
    colour = "black",
    linewidth = 0.2) +
  scale_fill_gradientn(
    colours = vector.map.palette,
    trans = "log1p",
    guide = "none") +
  scale_x_continuous(
    trans = "log1p",
    breaks = vector.map.breaks,
    limits = range(vector.map.breaks)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank())

# Add inset histogram
pred.map.with.hist <-
  pred.map +
  inset_element(
    hist.inset,
    left   = 0.45,
    bottom = 0.15,
    right  = 0.85,
    top    = 0.45)
pred.map.with.hist
# Save predicted vector count map as a pdf and as and RDS R object
ggsave(paste0("figures_public/", genus, "_pred_count_mean.pdf"), pred.map.with.hist, 
       height = 4, width = 7)
saveRDS(pred.map.with.hist, paste0("figures_public/", genus, "_pred_count_mean.rds"))

# ### Prepare publication figures

# #### Composite figure generation

# Prepare and export composite figures for publication
forest.best <- forest.plots[[Best]] + labs(tag = LETTERS[7 + (genus == "Aedes")])
forest.best@labels$title <- ""
forest.best

# Save forest plot as RDS file to be combined later
saveRDS(forest.best, paste0("figures_public/", genus, "_forest_plot.rds"))

# Load Culex and Aedes maps and forest plots and combine in 
# a single composite figure
vector.map.list.culex <- readRDS("figures_public/Culex_obs_count_mean.rds")
vector.map.list.aedes <- readRDS("figures_public/Aedes_obs_count_mean.rds")
forest.best.culex <- readRDS("figures_public/Culex_forest_plot.rds")
forest.best.aedes <- readRDS("figures_public/Aedes_forest_plot.rds")
count.map.plus.forest <-   
  wrap_plots(c(vector.map.list.culex, 
               forest.best.culex, 
               vector.map.list.aedes,
               forest.best.aedes), 
             byrow = FALSE, ncol = 2) +
    plot_layout(guides = "collect")

ggsave("figures_public/count_map_plus_forest.pdf", count.map.plus.forest, 
       height = 9, width = 11)

# #### Export prediction datasets

# Prepare the data set for export, including the mean and SD predictions
keep.cols <- 
  c("newID", "exposed", "District", "Village", "VillageTV", "hhID",
    "samplingRound", "ndvi.sds", "ndvi.sds2", "Class", 
    "NRainDays.log10.sds", "elevation500m.sds", 
    "buildingDensity500m.log10.sds", "popDens500m.log10.sds", 
    "proportionIrrigatedCrop.sds",
    paste0("mean", genus),
    "geometry")


# Convert geometry to WKT (well known text) format
pred.wkt <- mos.all[mos.all$Sampling == "Cattle", keep.cols]
pred.wkt$geometry <- st_as_text(st_geometry(pred.wkt))

# Step 2: Drop sf class (optional but cleaner)
pred.wkt <- as.data.frame(pred.wkt)

# Step 3: Write to CSV
write.csv(pred.wkt, 
          paste0("data_public/cattleSerologyWith", genus, "Pred.csv"), 
          row.names = FALSE)

# Stop timer
global.finish <- Sys.time()
print(global.finish - global.start)

