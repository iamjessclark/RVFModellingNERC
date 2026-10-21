#' Function to calculate Bayes factors from marginal likelihoods,
#' relative to the maximum marginal likelihood
bayes.factor <- function(mlik) exp(max(mlik) - mlik)

#' Function to extract DIC from a list of INLA fits, and calculate delta DIC
IC.inla <- 
  function(fit.list, crit = c("dic", "waic", "mlik")[1], Round = TRUE) {
    if(crit == "mlik") {
      tab <- 
        cbind.data.frame(IC = sapply(fit.list, function(f) 
          as.vector(f[[crit]]["log marginal-likelihood (integration)", ])))
      tab$BayesFactor <- bayes.factor(tab$IC)
    } else {
      tab <- cbind.data.frame(IC = sapply(fit.list, function(f) f[[crit]][[crit]]))
      tab$deltaIC <- tab$IC - min(tab$IC)
      tab$Weights <- exp(-0.5 * tab$deltaIC) / sum(exp(-0.5 * tab$deltaIC))
    }
    names(tab) <- gsub("IC", toupper(crit), names(tab))
    if(Round) tab <- round(tab, 2)
    tab
  }

# IC.inla(fit.list = Model.fits, crit = "mlik")
# IC.inla(fit.list = Model.fits, crit = "dic")
# IC.inla(fit.list = Model.fits, crit = "waic")

#' Function to fit a NA-response model to get quick approximate 
#' priors from inlabru and INLA fits
inla.na.fit <- 
  function(object, response.name = "BRU.response") {
    stopifnot(!is.null(object$.args$data[[response.name]]))
    object$.args$data[[response.name]] <- 
      rep(NA, length(object$.args$data[[response.name]]))
    inla.rerun(inla.rerun(object, plain = TRUE))
  }


#' Function to report messages (e.g. progress) from inside mclapply
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}

#' Slightly modified version of Greg Albery's function, with 
#' the discrete legend labels showing a range rather than the lower
#' bound
ggField <- 
  function (Model, Mesh, Groups = 1, GroupVar = NULL, GroupLabels = NULL, 
            Fill = "Discrete", NFillCat = 9,
            FillAlpha = F, FitScale = "Link", 
            LegendScaleFn = I, GrandMean = NULL, 
            GrandSD = NULL, Round = "Before", Boundary = NULL, BoundaryWidth = 1, 
            Res = 300, Points = NULL, PointAlpha = 1, PointSize = NULL, 
            PointColour = "black") {
    require(ggplot2)
    require(INLA)
    require(tidyverse)
    Projection <- inla.mesh.projector(Mesh, dims = c(Res, Res))
    Full.Projection <- expand.grid(x = Projection$x, y = Projection$y)
    Dim1 <- nrow(Full.Projection)
    WName <- Model$summary.hyperpar %>% rownames
    WName <- WName[str_detect(WName, "Range for")] %>% str_split(" ") %>% 
      map_chr(last)
    if (Groups == 1) {
      Full.Projection$value <- c(inla.mesh.project(Projection, 
                                                   Model$summary.random$w$mean))
    } else {
      Full.Projection[, paste0("Group", 1:Groups)] <- apply(matrix(Model$summary.random[[WName]]$mean, 
                                                                   ncol = Groups), 2, function(x) c(inla.mesh.project(Projection, 
                                                                                                                      x)))
      Full.Projection <- reshape2::melt(Full.Projection, id.vars = c(names(Full.Projection)[-which(names(Full.Projection) %in% 
                                                                                                     paste0("Group", 1:Groups))]))
    }
    Full.Projection$Group <- rep(1:Groups, each = Dim1)
    if (length(GrandMean) == Groups) {
      GrandMean %<>% rep(each = Dim1)
    }
    if (FitScale == "Binomial" & Round == "Before") {
      Full.Projection %<>% mutate_at("value", ~logistic(as.numeric(as.character(.x)) + 
                                                          GrandMean))
    }
    if (FitScale == "Gaussian" & Round == "Before") {
      Full.Projection %<>% mutate_at("value", ~(as.numeric(as.character(.x))) * 
                                       GrandSD + GrandMean)
    }
    if (Fill == "Discrete") {
      field.breaks <- quantile(Full.Projection$value, 0:NFillCat/NFillCat, na.rm = T)
      ndp <- ceiling(-log10(diff(range(LegendScaleFn(field.breaks)/100))))
      field.breaks.labels <-
        paste(round(LegendScaleFn(field.breaks), ndp)[1:NFillCat], 
              round(LegendScaleFn(field.breaks), ndp)[1:NFillCat + 1], sep = ", ")
      Full.Projection$Fill <-
        cut(Full.Projection$value,
            breaks = field.breaks, 
            labels = field.breaks.labels, 
            include.lowest = T)
    }
    else {
      Full.Projection$Fill <- Full.Projection$value
    }
    if (FitScale == "Gaussian" & Round == "After") {
      Full.Projection %<>% mutate_at("Fill", ~(as.numeric(as.character(.x))) * 
                                       GrandSD + GrandMean)
    }
    if (FitScale == "Binomial" & Round == "After") {
      Full.Projection %<>% mutate_at("value", ~logistic(as.numeric(as.character(.x)) + 
                                                          GrandMean))
      if (Fill == "Discrete") {
        Full.Projection %<>% mutate_at("Fill", ~factor(.x, 
                                                       levels = gtools::mixedsort(unique(.x))))
      }
    }
    Full.Projection <- na.omit(Full.Projection)
    FieldPlot <- ggplot(Full.Projection, aes(x, y))
    if (!is.null(Boundary)) {
      names(Boundary)[1] <- "x"
      names(Boundary)[2] <- "y"
      FieldPlot <- FieldPlot + geom_polygon(data = Boundary, 
                                            fill = "white")
    }
    if (!FillAlpha) {
      FieldPlot <- FieldPlot + geom_tile(colour = NA, aes(fill = Fill)) + 
        coord_fixed() + labs(fill = "Mean") + labs(x = "Easting", 
                                                   y = "Northing")
    }
    else {
      FieldPlot <- FieldPlot + geom_tile(colour = NA, aes(fill = Fill, 
                                                          alpha = Fill)) + coord_fixed() + labs(fill = "Mean") + 
        labs(x = "Easting", y = "Northing")
    }
    if (Fill == "Discrete") {
      FieldPlot <- FieldPlot + guides(fill = guide_legend(reverse = T))
    }
    if (!is.null(Boundary)) {
      FieldPlot <- FieldPlot + geom_polygon(data = Boundary, 
                                            fill = NA, colour = "black", size = BoundaryWidth, 
                                            inherit.aes = F, aes(x, y))
    }
    if (Groups > 1) {
      if (!is.null(GroupLabels)) {
        FacetLabels <- GroupLabels
      }
      else FacetLabels <- as.character(1:Groups)
      names(FacetLabels) <- as.character(1:Groups)
      FieldPlot <- FieldPlot + facet_wrap(~Group, labeller = as_labeller(FacetLabels)) + 
        theme(strip.background = element_rect(fill = "white"))
    }
    if (!is.null(Points)) {
      Points$X <- Points[, 1]
      Points$Y <- Points[, 2]
      if (Groups == 1) {
        FieldPlot <- FieldPlot + geom_point(data = Points, 
                                            inherit.aes = F, aes(X, Y), alpha = PointAlpha, 
                                            colour = PointColour)
      }
      else {
        Points$Group <- Points[, GroupVar] %>% as.factor %>% 
          as.numeric
        FieldPlot <- FieldPlot + geom_point(data = Points, 
                                            inherit.aes = F, aes(X, Y, group = Group), alpha = PointAlpha, 
                                            colour = PointColour, size = PointSize)
      }
    }
    return(FieldPlot)
  }


#' See example at inla.doc("gpoisson")
dgpoisson <- function(y, mu, phi, p) {
  a = mu + phi * mu^(p-1.0) * y;
  b = 1. + phi * mu^(p-1.0);
  d = exp(log(mu) + (y-1.0)*log(a) -
            y*log(b) - lfactorial(y) - a/b)
  return (d)
  2
}

#' See example at inla.doc("gpoisson")
rgpoisson <-
  function(n, mu, phi, p) {
  stopifnot(length(mu) == 1)
  s = sqrt(mu*(1+phi*mu^(p-1))^2)
  f = 20
  low = as.integer(max(0, mu - f*s))
  high = as.integer(mu + f*s)
  prob = dgpoisson(low:high, mu, phi, p)
  y = sample(low:high, n, replace=TRUE,
             prob = prob)
  return (y)
}


variance.decomp <-
  function(result, mm, spatial.re, random.effects = NULL, fixed.est.mean = NULL,
           re.var.stat = "mean", nb.ol.var.method = c("delta", "trigamma", "exclude")[1]) {
    
    distribution <- result$.args$family
    
    #' get variance components 
    #' fixed effects variance
    if(is.null(fixed.est.mean)) fixed.est.mean <- result$summary.fixed$mean
    if(length(mm) == 1 && is.na(mm)) {
      mu <- fixed.est.mean
      fixed.var <- 0
    } else {
      fixed.predictions <- as.matrix(mm) %*% cbind(fixed.est.mean)
      mu <- exp(fixed.predictions)
      fixed.var <- var(fixed.predictions)
    }
    #' IID random effects variance
    if(!is.null(random.effects)) {
      re.var <-
        sapply(random.effects, function(re) {
          1/result$summary.hyperpar[paste("Precision for", re), re.var.stat]
        })
    } else re.var <- NA
    #' spatial random effect variance
    spatial.var <-
      result$summary.hyperpar[paste("Stdev for", spatial.re), re.var.stat]^2
    #' observation-level variance on log scale 
    #' (see https://doi.org/10.1098/rsif.2017.0213)
    if(distribution == "nbinomial") {
      nb.theta <-
        result$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", re.var.stat]
      if(nb.ol.var.method == "trigamma") {
        ol.var <- mean(trigamma(1/(1/mu + 1/nb.theta)))
      }
      if(nb.ol.var.method == "delta") {
        ol.var <- mean(1/mu + 1/nb.theta)
      }
      if(nb.ol.var.method == "exclude") {
        ol.var <- 0
        warning("Observation-level variance excluded (set to zero)")
      }
    }
    #' using delta method to derive observation-level variance on log scale 
    if(distribution == "gpoisson") {
      phi <-
        result$summary.hyperpar["Overdispersion for gpoisson", re.var.stat]
      p <-
        ifelse(result$all.hyper$family[[1]]$hyper$theta2$fixed, 
               result$all.hyper$family[[1]]$hyper$theta2$initial, 
               result$summary.hyperpar["Parameter p_intern for gpoisson", "mean"]) 
      ol.var <-
        mean(1/mu + 2 * phi * mu^(p - 2) + phi^2 * mu^(2 * p - 3))
      if(nb.ol.var.method == "exclude") {
        ol.var <- 0
        warning("Observation-level variance excluded (set to zero)")
      }
    }
    #' Tabulate variance components and calculate marginal and conditional R^2
    var.components <-
      cbind(c(Fixed_effects_var = fixed.var,
              RE_var.Spatial = spatial.var,
              RE_var = re.var,
              Observation_var = ol.var))
    R2.marginal <- as.vector(var.components["Fixed_effects_var", ]/sum(var.components, na.rm = TRUE))
    R2.conditional <- as.vector(1 - var.components["Observation_var", ] / sum(var.components, na.rm = TRUE))
    
    var.components <- rbind(var.components, 
                            Total_var = sum(var.components, na.rm = TRUE),
                            R2.marginal = R2.marginal, 
                            R2.conditional = R2.conditional)
    cbind(na.omit(var.components))
  }

# variance.decomp(result = Model.fits$Environmental, mm = X.list$Environmental, 
#                 spatial.re = "field.tve",
#                 fixed.est.mean = Model.fits$Environmental$summary.random$X$mean,
#                 random.effects = c("Village", "VillageTV"),
#                 nb.ol.var.method = "delta")


# n <- 100000
# mu <- 50
# phi <- 0.056
# p <- 2.5
# nb.theta <- 1
# 
# x <- rnbinom(n, size = nb.theta, mu = mu)
# var(x)
# mu + mu^2 / nb.theta
# 
# trigamma(1/(1/mu + 1/nb.theta))
# 1/mu + 1/nb.theta
# var(log(x + 0.1))
# 
# y <- rgpoisson(n, mu = mu, phi = phi, p = p)
# hist(y)
# var(y)
# mu * (1 + 2 * phi * mu^(p - 1) + phi^2 * mu^(2 * p - 2))
# mu + 2 * phi * mu^p + phi^2 * mu^(2 * p - 1)
# 1/mu + 2 * phi * mu^(p - 2) + phi^2 * mu^(2 * p - 3)
# var(log(y + 0.1))


clean_one <- function(g) {
  r23 <- g[!is.na(g$status2023), ]
  r24 <- g[!is.na(g$status2024), ]
  
  # --- clean 2023 + 2024 pair ---
  if (nrow(g) == 2 && nrow(r23) == 1 && nrow(r24) == 1) {
    age_diff <- abs(r23$newAgeInt - r24$newAgeInt)
    if (is.na(age_diff)) {
      g$resolution <- "review: newAgeInt is NA"
      return(g)
    }
    if (age_diff > 1) return(g[0, ])                 # too far apart -> drop both
    merged <- r23
    for (col in names(g)) {
      if (!is.list(g[[col]])) merged[[col]] <- coalesce(r23[[col]], r24[[col]])
    }
    merged$newAgeInt  <- r23$newAgeInt
    merged$changes    <- "merged 2023+2024"
    merged$resolution <- paste0("merged (age diff ", age_diff, ")")
    return(merged)
  }
  
  # --- NEW: no 2023 row, every row is a 2024 entry ---
  if (nrow(r23) == 0 && nrow(r24) == nrow(g)) {
    if (n_distinct(r24$status2024) == 1) {           # status2024 agrees
      out <- g[1, ]                                  #   -> keep first occurrence
      out$resolution <- "kept first (matching status2024)"
      return(out)
    } else {                                         # status2024 differs
      return(g[0, ])                                 #   -> drop both
    }
  }
  
  # --- anything else still unexpected ---
  g$resolution <- "review: not a clean 2023/2024 pair"
  return(g)
}

