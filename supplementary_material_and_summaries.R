#' Summary and supplementary figures. Run FOI_main_analysis.R first in the same session.

library(tidytext)
library(gridExtra)

path = "FOI_analysis/public_outputs/"

################################################################################
#### Data exploration ####
################################################################################

hist(dat$meanCulex)
hist(dat$culexAbundance.sds)

hist(dat$meanAedes)
hist(dat$aedesAbundance.sds)

# how correlated is everything?
plot(st_drop_geometry(dat[, c("buildingDensity500m.log10.sds", "popDens500m.log10.sds", "elevation500m.sds", "proportionIrrigatedCrop.sds", "NRainDays.log10.sds", "ndvi.sds", "culexAbundance.sds", "aedesAbundance.sds")]))

################################################################################
#### Sampling summaries ####
################################################################################

cattle.age %>%
  group_by(District, hhID) %>%
  summarise(n = length(cattle_ID), .groups = "drop") %>%
  ggplot() +
  geom_col(aes(x = reorder_within(hhID, n, District), y = n)) +
  scale_x_reordered() +
  facet_grid(. ~ District, scales = "free_x")+
  xlab("ignore - just household IDs") +
  ylab("n cattle per household")


cattle.age %>%
  group_by(District, Village, classification) %>%
  summarise(n = length(unique(hhID)), .groups = "drop") %>%
  ggplot() +
  geom_col(aes(x = reorder_within(Village, n, District), y = n, fill = classification)) +
  scale_x_reordered() +
  facet_grid(. ~ District, scales = "free_x")+
  xlab("village") +
  ylab("n households per village")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


cattle.age %>%
  group_by(District, Village, classification) %>%
  summarise(n = length(unique(cattle_ID)), .groups = "drop") %>%
  ggplot() +
  geom_col(aes(x = reorder_within(Village, n, District), y = n, fill = classification)) +
  scale_x_reordered() +
  facet_grid(. ~ District, scales = "free_x")+
  xlab("village") +
  ylab("n cattle per village")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

cattle.age %>%
  group_by(newAgeInt) %>%
  tally() %>%
  filter(!is.na(newAgeInt)) %>%
  ggplot() +
  geom_col(aes(x = as.factor(newAgeInt), y = n)) +
  xlab("age") +
  theme_bw()
ggsave(path = path, "agePlot.pdf")

cattle.age %>%
  group_by(District, Village, newAgeInt) %>%
  tally() %>%
  ggplot() +
  geom_col(aes(x = Village, y = n, fill = newAgeInt)) +
  xlab("Village ID") +
  theme_bw() +
  coord_flip() +
  facet_grid(District ~., scales = "free") +
  scale_fill_viridis_c(breaks = seq(from = 1, to = 15, by = 1))+
  labs(fill = "Age (years)") +
  guides(fill = guide_legend(
    title.position = "top",
    title.hjust = 0.5
  ))
ggsave(path = path, "villageN.pdf")

################################################################################
#### Seroprevalence ####
################################################################################

##### Seroprevalence by age, district and year #####

cattle.age %>%
  dplyr::select(District, Village, newAgeInt, status2023, status2024) %>%

  group_by(District) %>%
  mutate(villNum = dense_rank(Village)) %>%
  ungroup() %>%

  pivot_longer(
    cols = c(status2023, status2024),
    names_to = "Year",
    values_to = "status"
  ) %>%

  filter(!is.na(status)) %>%

  count(District, newAgeInt, Year, status, name = "n") %>%

  group_by(District, newAgeInt, Year) %>%
  mutate(prevalence = n / sum(n) * 100,
         Year = as.factor(Year)) %>%
  ungroup() %>%

  mutate(Year = recode(Year,
                       status2023 = "2023",
                       status2024 = "2024")) %>%
  ggplot(aes(x = factor(newAgeInt), y = prevalence, fill = status)) +
  geom_col(position = "dodge", alpha = 0.6) +
  geom_text(aes(label = paste0(round(prevalence,0), "%", " (", n, ")")),
            position = position_dodge(width = 1),
            hjust = -1, size = 2.5, show.legend = FALSE) +
  xlab("Age (years)") +
  theme_bw() +
  scale_fill_viridis_d() +
  coord_flip() +
  facet_grid(District ~ Year, scales = "free")
ggsave(path = path, "expStatusNY.pdf")

##### Seroprevalence by year #####

cattle.age %>%
  dplyr::select(District, Village, newAgeInt, status2023, status2024) %>%

  group_by(District) %>%
  mutate(villNum = dense_rank(Village)) %>%
  ungroup() %>%

  pivot_longer(
    cols = c(status2023, status2024),
    names_to = "Year",
    values_to = "status"
  ) %>%

  filter(!is.na(status)) %>%

  count(Year, status, name = "n") %>%

  group_by(Year) %>%
  mutate(prevalence = n / sum(n) * 100) %>%
  ungroup() %>%

  mutate(Year = recode(Year,
                       status2023 = "2023",
                       status2024 = "2024")) %>%
  ggplot(aes(x = Year, y = prevalence, fill = status)) +
  geom_col(position = "dodge", alpha = 0.6) +
  geom_text(aes(label = paste0(round(prevalence,0), "%", " (", n, ")")),
            position = position_dodge(width = 1),
            vjust = -1, size = 2.5, show.legend = FALSE) +
  xlab("Age (years)") +
  theme_bw() +
  scale_fill_viridis_d()
ggsave(path = path, "expStatusNYear.pdf")

##### Seroprevalence by age within each village #####

cattle.age$exposed2023 <- ifelse(cattle.age$status2023 == "Positive", 1, 0)
cattle.age$exposed2024 <- ifelse(cattle.age$status2024 == "Positive", 1, 0)

cattle.age %>%
  dplyr::select(Village, newAgeInt, exposed2023, exposed2024) %>%
  pivot_longer(
    cols      = c(exposed2023, exposed2024),
    names_to  = "year",
    values_to = "status"
  ) %>%
  group_by(Village, newAgeInt) %>%
  summarise(
    n          = n(),
    cases      = sum(status, na.rm=T),
    prevalence = mean(status, na.rm=T),          # proportion of 1s
    .groups    = "drop"
  ) %>%
  ggplot() +
  geom_point(aes(x = newAgeInt, y = prevalence*100)) +
  facet_wrap(.~Village, scales = "free")

##### Village seroprevalence maps (2023 and 2024) #####

#' Give each District a core colour. Reused later for the facet strip
#' colours on the village-seroprevalence-by-district-and-year figure.
districtCoreColours <-
  setNames(RColorBrewer::brewer.pal(nlevels(dat$District), "Dark2"),
           levels(dat$District))

#' Build a village-level seroprevalence map (hue = District, tone = %) plus
#' its matching legend, for any binary 0/1 response column in dat. Refactored
#' into a function so the 2023 and 2024 maps share identical logic rather
#' than risking two hand-edited copies drifting out of sync.
buildSeroMap <- function(dat, responseVar, districtCoreColours, caption,
                          outlierThreshold = Inf) {

  #' Village-level seroprevalence (% positive), used to set the *tone* within
  #' each district's core colour - i.e. hue = District, tone = seroprevalence.
  #' n counts only non-NA responses (exposed24 is NA for animals already
  #' positive in 2023); villages with no data for this response are dropped.
  villageSeroprev <-
    dat %>%
    st_drop_geometry() %>%
    group_by(District, Village) %>%
    summarise(seroprevalence = mean(.data[[responseVar]], na.rm = TRUE) * 100,
              n = sum(!is.na(.data[[responseVar]])),
              .groups = "drop") %>%
    filter(n > 0)

  #' Villages at/above outlierThreshold (e.g. a single 100% village off a
  #' tiny sample) are excluded from the colour scale and shown as solid
  #' black instead, so they don't compress every other village's colour
  #' into an indistinguishable pale cluster
  isOutlier <- villageSeroprev$seroprevalence >= outlierThreshold

  #' Cap the colour scale at the highest *non-outlier* seroprevalence
  #' observed (rather than a fixed 0-100%), so the full colour range is used
  maxSeroprev <- max(villageSeroprev$seroprevalence[!isOutlier], na.rm = TRUE)

  #' Blend each district's core colour from pale (0% seroprevalence) to the
  #' full core colour (at maxSeroprev), then look up each village's tone;
  #' outlier villages get solid black regardless of district
  villageSeroprev$villageColour <-
    vapply(seq_len(nrow(villageSeroprev)), function(i) {
      if (isOutlier[i]) return("black")
      rampFun <- colorRamp(c("grey90", districtCoreColours[[villageSeroprev$District[i]]]))
      rgbVals <- rampFun(pmin(villageSeroprev$seroprevalence[i] / maxSeroprev, 1))
      rgb(rgbVals[1], rgbVals[2], rgbVals[3], maxColorValue = 255)
    }, character(1))

  #' One point per village (its centroid), rather than one per animal/
  #' household. inner_join drops villages with no data for this response.
  villageMap <-
    dat %>%
    group_by(District, Village) %>%
    summarise(.groups = "drop") %>%
    st_centroid() %>%
    inner_join(villageSeroprev, by = c("District", "Village"))

  #' Small legend showing the seroprevalence-to-tone mapping for each district,
  #' scaled 0 to maxSeroprev (matching the colours actually used above).
  #' Discrete swatches every 10 percentage points (plus the exact max), with
  #' gaps between them, rather than one smooth blended bar - much easier to
  #' read off a specific value when the max is a less-than-100 number like 52%.
  #' Drop a round-10 break if it's within 3 points of the max, to avoid the
  #' two labels (e.g. "50%" and "52%") crowding into each other
  seroTenBreaks <- seq(0, floor(maxSeroprev / 10) * 10, by = 10)
  seroTenBreaks <- seroTenBreaks[abs(seroTenBreaks - maxSeroprev) > 3]
  seroBreaks <- c(seroTenBreaks, round(maxSeroprev))
  #' Size each swatch to *its own* local spacing (not one global width) - the
  #' regular 10-point gaps get nice fat swatches, and only the short final
  #' segment (round(maxSeroprev) sitting close to the last multiple of 10)
  #' gets a correspondingly thinner one, rather than shrinking everything.
  seroGapBefore <- c(Inf, diff(seroBreaks))
  seroGapAfter <- c(diff(seroBreaks), Inf)
  seroTileWidths <- pmin(seroGapBefore, seroGapAfter) * 0.7

  seroLegend <-
    lapply(names(districtCoreColours), function(d) {
      rampFun <- colorRamp(c("grey90", districtCoreColours[[d]]))
      #' pmin(..., 1): round(maxSeroprev) can round up past the true
      #' maxSeroprev (e.g. 52.6 -> 53), which would push the ratio above 1 -
      #' colorRamp() returns NA (not a clamped colour) for anything outside
      #' [0,1], which is what was throwing the "color intensity NA" error
      data.frame(x = seroBreaks,
                 width = seroTileWidths,
                 colour = rgb(rampFun(pmin(seroBreaks / maxSeroprev, 1)), maxColorValue = 255)) %>%
        ggplot(aes(x = x, y = 1, fill = colour, width = width)) +
        geom_tile() +
        scale_fill_identity() +
        scale_x_continuous(breaks = seroBreaks, labels = paste0(seroBreaks, "%"),
                           expand = expansion(add = max(seroTileWidths))) +
        ggtitle(d) +
        theme_void(base_size = 8) +
        theme(axis.text.x = element_text(size = 6, vjust = 3),
              plot.title = element_text(size = 8, hjust = 0.5, colour = districtCoreColours[[d]], face = "bold"))
    })

  #' If any village was excluded from the scale as an outlier, add a small
  #' black swatch to the legend so it isn't just an unexplained black dot
  if (any(isOutlier)) {
    outlierLabel <- if (min(villageSeroprev$seroprevalence[isOutlier]) >= 100) {
      "100%"
    } else {
      paste0("≥", round(min(villageSeroprev$seroprevalence[isOutlier])), "%")
    }
    outlierLegend <-
      data.frame(x = 1, colour = "black") %>%
      ggplot(aes(x = x, y = 1, fill = colour)) +
      geom_tile(width = 0.7) +
      scale_fill_identity() +
      scale_x_continuous(breaks = 1, labels = outlierLabel, limits = c(0.4, 1.6)) +
      ggtitle("Outlier") +
      theme_void(base_size = 8) +
      theme(axis.text.x = element_text(size = 6, vjust = 3),
            plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))
    seroLegend <- c(seroLegend, list(outlierLegend))
  }

  # Map of the sample sites, one point per village: hue = District, tone = seroprevalence
  #' shape 21/22/24 are the "fillable" point shapes (circle/square/triangle)
  #' that support colour (border) and fill (interior) separately - a plain
  #' shape=21 literal, or aes(colour=...), can't give a dark-grey border with
  #' a coloured interior at the same time, which is why fill needs its own aes
  siteMap <-
    ggplot(villageMap, aes(shape = District, fill = villageColour)) +
    annotation_map_tile(type = "cartolight", # available tile types: rosm::osm.types()
                        zoomin = 0) +
    geom_sf(inherit.aes = TRUE,
            size = 3,
            colour = "darkgrey") +
    scale_fill_identity() +
    scale_shape_manual(values = c(21, 22, 24)) +
    annotation_scale(location = "br", width_hint = 0.25, height = unit(0.1, "cm")) +
    theme_minimal() +
    labs(caption = caption)

  #' Use a relative height ratio, not an absolute unit for the legend row - the
  #' map has a locked geographic aspect ratio (via coord_sf/annotation_map_tile),
  #' and mixing that with a fixed "cm" row fights the lock and can squash the
  #' map down to almost nothing. A generous relative ratio keeps the legend
  #' visually short while staying robust to the map's own aspect constraint.
  combined <- siteMap / wrap_plots(seroLegend, ncol = length(seroLegend)) + plot_layout(heights = c(12, 1))

  list(villageMap = villageMap, siteMap = siteMap, seroLegend = seroLegend, combined = combined)
}

seroMap2023 <- buildSeroMap(dat, "exposed", districtCoreColours,
                            "village-level seroprevalence (2023)")
seroMap2023$combined

#' outlierThreshold = 99: the Babati village at 100% seroconversion is shown
#' in solid black rather than setting the colour scale, so the other, much
#' lower villages remain visually distinguishable from each other
seroMap2024 <- buildSeroMap(dat, "exposed24", districtCoreColours,
                            "village-level seroconversion (2024)",
                            outlierThreshold = 99)
seroMap2024$combined

#' Export village-level seroprevalence/seroconversion for both years as one
#' long table (one row per village per year). Villages with no 2024 data
#' (e.g. every animal already seropositive in 2023) simply have no 2024 row.
villageSeroprevByYear <-
  bind_rows(
    seroMap2023$villageMap %>% st_drop_geometry() %>% mutate(Year = 2023),
    seroMap2024$villageMap %>% st_drop_geometry() %>% mutate(Year = 2024)
  ) %>%
  mutate(seroprevalence = round(seroprevalence, 1)) %>%
  dplyr::select(District, Village, Year, seroprevalence, n) %>%
  arrange(District, Village, Year)

write.csv(villageSeroprevByYear, "FOI_analysis/public_outputs/villageSeroprevalenceByYear.csv", row.names = FALSE)

##### Village seroprevalence by district and year #####

villageSeroByYear <-
  cattle.age %>%
  dplyr::select(District, Village, status2023, status2024) %>%

  pivot_longer(
    cols = c(status2023, status2024),
    names_to = "Year",
    values_to = "status"
  ) %>%

  filter(!is.na(status)) %>%

  mutate(status = if_else(status == "Doubtful", "Negative", status)) %>%

  count(District, Village, Year, status, name = "n") %>%

  group_by(District, Village, Year) %>%
  mutate(prevalence = n / sum(n) * 100) %>%
  ungroup() %>%

  mutate(Year = recode(Year,
                       status2023 = "2023",
                       status2024 = "2024")) %>%

  filter(status == "Positive")

#' ggplot's facet_wrap() can't colour individual strips differently, so build
#' one facetted-by-Year plot per District (each with its own strip colour,
#' matching districtCoreColours from the site map above) and combine them
#' with patchwork into the same ncol = 3 layout as the original facet_wrap
districtPlots <-
  lapply(names(districtCoreColours), function(d) {
    villageSeroByYear %>%
      filter(District == d) %>%
      ggplot(aes(x = reorder_within(Village, -prevalence, Year), y = prevalence)) +
      geom_col(alpha = 0.6, fill = districtCoreColours[[d]]) +
      geom_text(aes(label = paste0(round(prevalence, 0), "% (", n, ")")),
                vjust = -0.5, size = 2.5) +
      scale_x_reordered() +
      facet_wrap(~ Year, ncol = 1, scales = "free_x") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.background = element_rect(fill = districtCoreColours[[d]]),
            strip.text = element_text(colour = "white", face = "bold")) +
      xlab("Village") +
      ggtitle(d)
  })

wrap_plots(districtPlots, ncol = 3)
ggsave(path = path, "expStatusNDistrictYear.pdf")

##### Village seroprevalence coloured by model-estimated FOI #####

dat %>%
  st_drop_geometry() %>%
  group_by(Village) %>%
  mutate(group_mean = mean(villageFOI, na.rm=T)) %>%
  ungroup() %>%
  group_by(Village, group_mean) %>%
  summarise(
    n          = n(),
    cases      = sum(exposed),
    prevalence = mean(exposed),          # proportion of 1s
    .groups    = "drop"
  ) %>%
  ggplot() +
  geom_col(aes(x = reorder(Village, -prevalence), y = prevalence, fill = group_mean*100))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

##### Map of observed seroconversions 2023-2024 #####

villageSeroconversions <-
  cattleDatInc %>%
  mutate(X = st_coordinates(geometry)[, 1],
         Y = st_coordinates(geometry)[, 2]) %>%
  st_drop_geometry() %>%
  group_by(Village) %>%
  summarise(X = mean(X), Y = mean(Y), n = sum(exposed == 1, na.rm = TRUE)) %>%
  st_as_sf(coords = c("X", "Y"), crs = st_crs(cattleDatInc), remove = FALSE)

ggplot(villageSeroconversions) +
  annotation_map_tile(type = "cartolight", zoomin = 0) +  # You can change the tile type
  geom_sf(aes(fill = n),
          size = 3, shape = 21, colour = "darkgrey") +
  scale_fill_viridis_c( option = "G",
    direction = -1, breaks = scales::breaks_width(2)) +
  annotation_scale(location = "tl", width_hint = 0.5) +
  coord_sf(crs = crs.km, expand = TRUE) +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "bottom"
  ) +
  labs(fill = "n seroconversions")

ggsave("FOI_analysis/public_outputs/nseroConversions.pdf")

################################################################################
#### Model diagnostics ####
################################################################################

#' Plot the mesh
mesh.sfc <- fm_as_sfc(mods$inputs$spatial$mesh)
st_crs(mesh.sfc) <- st_crs(dat)
ggplot(data = dat) +
  geom_sf(data = mesh.sfc) +
  geom_sf(aes(color = Village)) +
  theme_minimal() +
  theme(legend.position="none")

bru_convergence_plot(mods$outputs$spatialVectors)

##### Cross-validation folds #####

cv_plot(cv, dat)

#' Make a map of the sampling locations showing folds
ggplot() +
  annotation_map_tile(type = "cartolight", # available types: rosm::osm.types()
                      zoomin = 0) +  # You can change the tile type
  geom_sf(data = dat.all[!is.na(dat.all$fold), ],
          aes(colour = as.factor(fold)),
          inherit.aes = FALSE,
          alpha = 0.6,
          size = 2) +
  annotation_scale(location = "tl", width_hint = 0.25, height = unit(0.1, "cm")) +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(colour = "Fold")

ggsave("FOI_analysis/public_outputs/kfoldblocks.pdf")

##### Village and herd FOI #####

#' Plot Village and herd
plot(herdFOI ~ villageFOI, data = dat, log = "xy")
abline(0, 1)
boxplot(villageFOI ~ Village, data = dat, las = 3)
boxplot(herdFOI ~ Village, data = dat, las = 3)

##### Village FOI against vector abundance #####

plot(FOI ~ culexAbundance, data = VillagePred)
cor.test(VillagePred$obsPos, VillagePred$culexAbundance)

plot(FOI ~ aedesAbundance, data = VillagePred)
cor.test(VillagePred$obsPos, VillagePred$aedesAbundance)

##### Observed vs predicted seropositives per village #####

ggplot(data = VillagePred,
       mapping = aes(x = predPos, y = obsPos, colour = meanAge, size = NTested)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point() +
  ylab("N observed seropositive per Village") +
  xlab("N predicted seropositive per Village") +
  scale_colour_viridis_c()

##### Observed vs predicted seroprevalence by age #####
dat$newAgeIntCat <- factor(paste0(dat$newAgeInt, "y"))

AgePred <-
  data.frame(Age = sort(unique(dat$newAgeInt)))
AgePred$AgeCat <- paste0(AgePred$Age, "y")
AgePred$AgeCat <- factor(AgePred$AgeCat, AgePred$AgeCat)
rownames(AgePred) <- AgePred$AgeCat
AgePred$NTested <-
  as.vector(table(dat$newAgeIntCat)[as.character(AgePred$AgeCat)])
AgePred$obsPos <-
  tapply(dat$exposed, dat$newAgeIntCat, sum)[as.character(AgePred$AgeCat)]
AgePred$predPos <-
  tapply(dat$probPos, dat$newAgeIntCat, sum)[as.character(AgePred$AgeCat)]

AgePred$obsPrev <- AgePred$obsPos/AgePred$NTested
AgePred$predPrev <- AgePred$predPos/AgePred$NTested

ggplot(data = AgePred,
       mapping = aes(x = Age, y = predPrev, size = NTested, colour = "black")) +
  geom_point(alpha = 0.5) +
  geom_point(aes(y = obsPrev, colour = "red"), alpha = 0.5) +
  scale_x_continuous(breaks = 1:15) +
  scale_size_continuous(breaks = c(1, 5, 10, 50, 100, 500, 1000),
                        name = "N tested") +
  ylab("Seroprevalence") +
  xlab("Age (years)") +
  scale_colour_manual(name = '',
                      values =c(black='black',red='red'),
                      labels = c('Predicted','Observed')) +
  theme_minimal()

################################################################################
#### Supplementary analysis: variation in FOI by year (latent RW1) ####
################################################################################

# Define the range of years
head(dat)
# Define the birth year of 1-year-olds as year 0 and previous years as -1, -2, etc
dat$birthyear <- 1 - dat$newAgeInt
years <- min(dat$birthyear):0
nyears <- length(years)


# Build a sparse matrix that links the years exposed to
# each individual. For each row (i) the matrix will have the value x=1
# for the corresponding columns indexed by j.
# Create a list of indices for each individual
A.list <-
  lapply(1:nrow(dat), function(i) {
    exposedyears <- dat$birthyear[i]:0
    match(exposedyears, years)
  })
# A.list indexes the exposed years for each animal. E.g. animal 1 is
dat$newAgeInt[1]
# years old, so was born in year
dat$birthyear[1]
# and was exposed for the following years
years[A.list[[1]]]
A.year <-
  sparseMatrix(
    i = rep(1:nrow(dat), sapply(A.list, length)), # the index of each individual repeated for each year of exposure (length = sum(ages))
    j = unlist(A.list), # the index of the years exposed for each individual (length = sum(ages))
    x = 1,
    dims = c(nrow(dat), nyears))


spatialA <- inla.spde.make.A(mods$inputs$spatial$mesh, loc = dat$geometry)
field <- inla.spde.make.index('field', n.spde =  mods$inputs$spatial$spde$n.spde)


StackHost.tv <- inla.stack(
  data = list(y = dat$exposed), # specify the response variable
  A = list(1, 1, 1, 1, spatialA, A.year), # Vector of Multiplication factors for random and fixed effects
  effects = list(
    Intercept = X[, 1],
    X = data.frame(newAgeInt = X[, 2]), # attach the model matrix
    Village = dat$Village,
    hhID = dat$hhID, # insert vectors of any random effects
    field = field,
    year.index = 1:nyears))

# From A.year, make one matrix for each district/cluster, and
# make a new StackYear for each.
dat$District <- factor(dat$District)
A.year.dist <-
  lapply(levels(dat$District), function(x) {
    Diagonal(x = as.numeric(dat$District == x)) %*% A.year
  })
names(A.year.dist) <- levels(dat$District)
lapply(A.year.dist, dim)

tapply(dat$newAgeInt, dat$District, max)

# add the year random effect to the model formula
model4 <- as.formula("y ~ -1 + Intercept + offset(log(newAgeInt)) +
                     f(field, model = mods$inputs$spatial$spde) +
                     f(year.index, model = 'rw1', hyper = mods$priors$iid.prec.prior) +
                     f(Village, model = 'iid', hyper = mods$priors$iid.prec.prior) +
                     f(hhID, model = 'iid', hyper = mods$priors$iid.prec.prior)")

# fit the model with a single global between-year random effect
fit4 <- inla(model4,
             family = "binomial",
             Ntrials = 1,
             control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
             control.family = list(link = "cloglog"),
             control.predictor = list(A = inla.stack.A(StackHost.tv), compute = TRUE),
             control.fixed = mods$priors$prior.fixed,
             data = inla.stack.data(StackHost.tv))
# the warning is fine - we want the identity link to be used

summary(fit4)

# ...the precision for the year random effect is very high, so there is basically
# no evidence for variation in FOI over years

ggplot(fit4$summary.random$Village, aes(y = ID, x = exp(mean))) +
  geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_y_discrete(limits = fit4$summary.random$Village$ID[order(fit4$summary.random$Village$mean)]) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
                     transform = "log") +
  theme_minimal() +
  labs(y = "", x = "Estimate",
       title = "Village effect estimates ± 95% CI")

X.dist <- as.data.frame(model.matrix(~ newAgeInt + District, data = dat))

StackHost.tv.dist <- inla.stack(
  data = list(y = dat$exposed), # specify the response variable
  A = list(1, 1, 1, 1, spatialA, A.year.dist$Babati, A.year.dist$Hai, A.year.dist$Moshi), # Vector of Multiplication factors for random and fixed effects
  effects = list(
    Intercept = X.dist[, 1],
    X = X.dist[, -1], # attach the model matrix
    Village = dat$Village,
    hhID = dat$hhID, # insert vectors of any random effects
    field = field,
    year.indexBabati = 1:nyears,
    year.indexHai = 1:nyears,
    year.indexMoshi = 1:nyears))



# add the district-specific year random effects to the model formula
model5 <- as.formula("y ~ -1 + Intercept + offset(log(newAgeInt)) +
                     DistrictHai + DistrictMoshi +
                     f(field, model = mods$inputs$spatial$spde) +
                     f(year.indexBabati, model = 'rw1', hyper = mods$priors$iid.prec.prior) +
                     f(year.indexHai, model = 'rw1', hyper = mods$priors$iid.prec.prior) +
                     f(year.indexMoshi, model = 'rw1', hyper = mods$priors$iid.prec.prior) +
                     f(Village, model = 'iid', hyper = mods$priors$iid.prec.prior) +
                     f(hhID, model = 'iid', hyper = mods$priors$iid.prec.prior)")


# fit the model with a single global between-year random effect
fit5 <- inla(model5,
             family = "binomial",
             Ntrials = 1,
             control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
             control.family = list(link = "cloglog"),
             control.predictor = list(A = inla.stack.A(StackHost.tv.dist), compute = TRUE),
             control.fixed = mods$priors$prior.fixed,
             data = inla.stack.data(StackHost.tv.dist))

# the warning is fine - we want the identity link to be used
summary(fit5)
GLMMmisc::mrr(1/fit5$summary.hyperpar["Precision for year.indexBabati", c("0.975quant", "0.025quant")])
GLMMmisc::mrr(1/fit5$summary.hyperpar["Precision for year.indexHai", c("0.975quant", "0.025quant")])
GLMMmisc::mrr(1/fit5$summary.hyperpar["Precision for year.indexMoshi", c("0.975quant", "0.025quant")])

# ...the precision for the year random effect is very high, so there is basically
# no evidence for variation in FOI over years

cowplot::plot_grid(

  ggplot(fit5$summary.random$year.indexBabati, aes(y = ID, x = exp(mean))) +
    geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_y_discrete(limits = fit5$summary.random$year.indexBabati$ID[order(fit5$summary.random$year.indexBabati$mean)]) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
                       transform = "log") +
    theme_minimal() +
    labs(y = "", x = "Estimate",
         title = "A"),

  ggplot(fit5$summary.random$year.indexHai, aes(y = ID, x = exp(mean))) +
    geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_y_discrete(limits = fit5$summary.random$year.indexHai$ID[order(fit5$summary.random$year.indexHai$mean)]) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
                       transform = "log") +
    theme_minimal() +
    labs(y = "", x = "Estimate",
         title = "B"),

  ggplot(fit5$summary.random$year.indexMoshi, aes(y = ID, x = exp(mean))) +
    geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_y_discrete(limits = fit5$summary.random$year.indexMoshi$ID[order(fit5$summary.random$year.indexMoshi$mean)]) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
                       transform = "log") +
    theme_minimal() +
    labs(y = "", x = "Estimate",
         title = "C"), ncol = 3, nrow = 1

)

# compare the models using DIC (dodgy) and WAIC (less dodgy)
modList <- list(SpatialYear = fit4, SpatialYearClust = fit5)
sapply(modList, function(f) f$waic$waic)
# ... the conclusion is the same: no evidence for a year effect on FOI

################################################################################
#### Vector counts by District and sampling round ####
################################################################################

vectorSummaryByDistrict <-
  read.csv("data_public/GenusLevel.csv") %>%
  select(District, samplingRound, Culex, Aedes) %>%
  group_by(District, samplingRound) %>%
  summarise(meanCulex = mean(Culex, na.rm = TRUE),
            sdCulex = sd(Culex, na.rm = TRUE),
            meanAedes = mean(Aedes, na.rm = TRUE),
            sdAedes = sd(Aedes, na.rm = TRUE),
            n = n(),
            .groups = "drop") %>%
  filter(District == "Babati" | District == "Hai" | District == "Moshi")

write.csv(vectorSummaryByDistrict, "FOI_analysis/public_outputs/vectorSummaryByDistrict.csv", row.names = FALSE)

#' Save the summary table as a pdf
vectorSummary.tab <-
  vectorSummaryByDistrict %>%
  mutate(across(c(meanCulex, sdCulex, meanAedes, sdAedes), ~round(.x, 2)))

pdf("FOI_analysis/public_outputs/vectorSummaryByDistrict.pdf", width = 8, height = 3)
grid.arrange(tableGrob(vectorSummary.tab, rows = NULL))
dev.off()
