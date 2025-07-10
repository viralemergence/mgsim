library(data.table)
library(here)
library(abc)
# library(triangle)
library(mc2d)
library(fitdistrplus)

options(scipen = 999)

demo_params <- fread(here("Data/Input/sample_data_round1.csv"))
abc_first_pass <- qread(here("Data/Validation/abc_round1.qs"))
round1_metrics <- read_csv(here("Data/Validation/round1_validation_metrics.csv")) |>
  mutate(hm_trend_1970 = rep(NA_real_), hm_trend_1993 = rep(NA_real_))

demo_params_selected <- demo_params[round1_metrics$dc,][, `:=`(Region = abc_first_pass$region,
                                           Weights = 1/(abc_first_pass$dist+.Machine$double.eps))][
                                             Region == TRUE, ]
summary(demo_params_selected)

demo_params_selected[1:10, ]

x <- demo_params_selected[["density_max"]]
w <- demo_params_selected[Region == "TRUE", ][["Weights"]]/sum(demo_params_selected[Region == "TRUE", ][["Weights"]])
min.max <- as.vector(quantile(x, c(0.05,0.95)))
{set.seed(54);density_x <- density(x, from = min.max[1], to = min.max[2],
                                   bw = "SJ", kernel = "gaussian",
                                   n = 1000, adjust = 0.25,
                                   weights = w)
  resample_x <- approx(
    x = cumsum(density_x$y)/sum(density_x$y),
    y = density_x$x,
    n = 10000,
    yleft = min(density_x$x), yright = max(density_x$x),
    na.rm = TRUE)$y
}
summary(resample_x)

par(mfrow = c(2, 1))
hist(x[x >= min.max[1] & x <= min.max[2]], breaks = 30)
hist(resample_x, breaks = 30)
par(mfrow = c(1,1))

# Truncated distribution functions

## normal
rnormt <- function(n, range, mean, sd, ...) {
  F.a <- pnorm(min(range), mean = mean, sd = sd, ...)
  F.b <- pnorm(max(range), mean = mean, sd = sd, ...)
  u <- runif(n, min = F.a, max = F.b)
  qnorm(p = u, mean = mean, sd = sd, ...)
}

## log-normal
rlnormt <- function(n, range, meanlog, sdlog, ...) {
  F.a <- plnorm(min(range), meanlog = meanlog, sdlog = sdlog, ...)
  F.b <- plnorm(max(range), meanlog = meanlog, sdlog = sdlog, ...)
  u <- runif(n, min = F.a, max = F.b)
  qlnorm(u, meanlog = meanlog, sdlog = sdlog, ...)
}

## gamma
rgammat <- function(n, range, shape, rate, ...) {
  F.a <- pgamma(min(range), shape = shape, rate = rate,...)
  F.b <- pgamma(max(range), shape = shape, rate = rate, ...)
  u <- runif(n, min = F.a, max = F.b)
  qgamma(u, shape = shape,rate = rate, ...)
}

## beta
rbetat <- function(n, range, shape1, shape2, ...) {
  F.a <- pbeta(min(range), shape1, shape2,...)
  F.b <- pbeta(max(range), shape1, shape2,...)
  u <- runif(n, min = F.a, max = F.b)
  qbeta(u, shape1, shape2,...)
}

## triangular
rtrianglet <- function(n, range, min, max, ...) {
  F.a <- ptriang(min(range), min = min, max = max, ...)
  F.b <- ptriang(max(range), min = min, max = max, ...)
  u <- runif(n, min = F.a, max = F.b)
  qtriang(u, min = min, max = max, ...)
}

## poisson
rpoist <- function(n, range, lambda, ...) {
  F.a <- ppois(min(range), lambda, ...)
  F.b <- ppois(max(range), lambda, ...)
  u <- runif(n, min = F.a, max = F.b)
  qpois(p = u, lambda, ...)
}

## negbin
rnbt <- function(n, range, size, mu, ...) {
  F.a <- pnbinom(min(range), size = size, mu = mu, ...)
  F.b <- pnbinom(max(range), size = size, mu = mu, ...)
  u <- runif(n, min = F.a, max = F.b)
  qnbinom(u, size = size, mu = mu, ...)
}

# Esimate parameters from the original 100 selected simulations
lnD <- fitdist(x, distr = "lnorm", method = "mge", gof = "KS")
triD <- fitdist(x, dist = "triang",
                start = list(min = min(min.max),
                             max = max(min.max),
                             mode = raster::modal(resample_x)),
                method = "mge", gof = "KS")
poisD <- fitdist(x, distr = "pois", method = "mme")
negbD <- fitdist(x, distr = "nbinom", method = "mme")

summary(lnD); summary(triD); summary(poisD); summary(negbD)
cdfcomp(list(lnD, triD, negbD))
qqcomp(list(lnD, triD, negbD))
ppcomp(list(lnD, triD, negbD))
gofstat(list(lnD, triD, negbD))

{plot(density(x, #from = min.max[1], to = min.max[2],
             bw = "SJ", kernel = "gaussian",
             n = 1000, adjust = 1,
             weights = w), ylim = c(0, 27))
abline(v = min.max, lty = 2, col = "grey80")
lines(density(resample_x, from = min(min.max), to = max(min.max)), lty = 2)
lines(density(rlnormt(n = 10000, range = min.max, meanlog = lnD$estimate[1], sdlog = lnD$estimate[2]),
              from = min(min.max), to = max(min.max)), col = "blue")
gRugs <- rgammat(n = 10000, range = min.max, shape = gD$estimate[1], rate = gD$estimate[2])
lines(density(gRugs,from = min(min.max), to = max(min.max)), col = "darkgreen")
bRugs <- rbetat(n = 10000, range = min.max, shape1 = bD$estimate[1], shape2 = bD$estimate[2])
lines(density(bRugs, from = min(min.max), to = max(min.max)), col = "red")
rug(x = bRugs, col = "red")
rug(x = x, col = "black")
rug(x = gRugs, side = 3, col = "darkgreen")
legend("topright",
       legend = c("Orig", "Resampled",
                  sprintf("lnorm = %s", round(lnD$aic)),
                  sprintf("gamma = %s", round(gD$aic)),
                  sprintf("beta = %s", round(bD$aic))),
       col = c(rep("black", 2), "blue", "darkgreen", "red"),
       lty = c(1, 2, rep(1, 3)),
       bty = "n")}

## Iterate through each demo variable and check distributions
param_dists <- lapply(colnames(demo_params_selected)[-c(25, 31, 48, 49)], function(param) {
  message(param)
  x <- demo_params_selected[[param]]
  w <- demo_params_selected[Region == "TRUE", ][["Weights"]]/sum(demo_params_selected[Region == "TRUE", ][["Weights"]])
  min.max <- as.vector(quantile(x, c(0.05,0.95)))
  # set seed and resample to 10,000 obs from ECDF
  ## https://stackoverflow.com/questions/32871602/r-generate-data-from-a-probability-density-distribution
  {set.seed(20212810);density_x <- density(x, from = min.max[1], to = min.max[2],
                                     bw = "SJ", kernel = "gaussian",
                                     n = 1000, adjust = 0.75,
                                     weights = w)
    resample_x <- approx(
      x = cumsum(density_x$y)/sum(density_x$y),
      y = density_x$x,
      n = 10000,
      yleft = min(density_x$x), yright = max(density_x$x),
      na.rm = TRUE)$y
  }
  par(mfrow = c(2, 1))
  {plot(density(x, from = min(x), to = max(x)), main = param)
     abline(v = min.max, col = "grey70", lty = 2, lwd = 2)}
  disc <- ifelse(param %in% c("density_max", "abundance_threshold"), TRUE, FALSE)
  if (param %in% c("density_max", "abundance_threshold")) {
    resample_x <- round(resample_x)
  }
  # Estimate parameters from interpolated/density estimated obs
  orig_dists <- list(
    normD <- tryCatch(fitdist(resample_x, distr = "norm", method = "mge",
                              gof = "KS", discrete = disc),
                      error = function(err) invisible({err})),
    lnD <- tryCatch(fitdist(resample_x, distr = "lnorm", method = "mge",
                            gof = "KS", discrete = disc),
                    error = function(err) invisible({err})),
    gD <- tryCatch(fitdist(resample_x, distr = "gamma", method = "mge",
                           gof = "KS", discrete = disc),
                   error = function(err) invisible({err})),
    bD <- tryCatch(fitdist(resample_x, distr = "beta", method = "mge",
                           gof = "KS", discrete = disc),
                   error = function(err) invisible({err})),
    triD <- tryCatch(fitdist(x, dist = "triang",
                             start = list(min = min(min.max),
                                          max = max(min.max),
                                          mode = raster::modal(x)),
                             method = "mge", gof = "KS", discrete = disc),
                     error = function(err) invisible({err})),
    unifD <- tryCatch(fitdist(resample_x, distr = "unif",method = "mge",
                              gof = "KS", discrete = disc),
                      error = function(err) invisible({err})),
    poisD <- tryCatch(fitdist(resample_x, distr = "pois", method = "mme",
                              discrete = TRUE),
                      error = function(err) invisible({err}),
                      warning = function(warn) invisible({warn})),
    negbD <- tryCatch(fitdist(resample_x, distr = "nbinom", method = "mme",
                              discrete = TRUE),
                      error = function(err) invisible({err}),
                      warning = function(warn) invisible({warn}))
  )
  # index for dist that errored
  keep_dists <- sapply(orig_dists, function(i) inherits(i, c("error", "warning")))
  orig_dists <- orig_dists[!keep_dists]
  names(orig_dists) <- c("norm", "lnorm", "gamma", "beta",
                         "triangular",
                         "unif", "pois", "negbin")[!keep_dists]
  cdfcomp(orig_dists, main = param, addlegend = TRUE, ylim = c(0, 1))
  par(mfrow = c(1,1))
  # gofstat(orig_dists, fitnames = names(orig_dists))
  return(
    list(
      Output = orig_dists)#,
      #GOF = gofstat(orig_dists, fitnames = names(orig_dists))
  )
})
names(param_dists) <- colnames(demo_params_selected)[-c(25, 31, 48, 49)]

# For each param, calculate the min AIC, then calc. delta AIC for each distribution
delta_aic <- lapply(names(param_dists), function(param) {
  var_dist <- param_dists[param][param][[1]]$Output
  dist_aic <- sapply(var_dist, function(i) i$aic)
  delta_aic <- dist_aic - min(dist_aic, na.rm = TRUE)
  delta_aic <- sort(delta_aic, na.last = TRUE)
  return(delta_aic)})
names(delta_aic) <- names(param_dists)
delta_aic

# Make 10,000 resamples based on results
{set.seed(20212511);
  dt_lhs_run3 <- data.table(
    ## standard deviation = beta
    standard_deviation = rbetat(n = 10000,
                                range = quantile(demo_params_selected[["standard_deviation"]], probs = c(0.05,0.95)),
                                shape1 = param_dists$standard_deviation$Output$beta$estimate[1],
                                shape2 = param_dists$standard_deviation$Output$beta$estimate[2]),
    ## growth rate max = log-norm
    growth_rate_max = rlnormt(n = 10000,
                              range = quantile(demo_params_selected[["growth_rate_max"]], probs = c(0.05,0.95)),
                              meanlog = param_dists$growth_rate_max$Output$lnorm$estimate[1],
                              sdlog = param_dists$growth_rate_max$Output$lnorm$estimate[2]),
    ## density max = normal
    density_max = round(
      rnormt(n = 10000,
             range = quantile(demo_params_selected[["density_max"]], probs = c(0.05,0.95)),
             mean = param_dists$density_max$Output$norm$estimate[1],
             sd = param_dists$density_max$Output$norm$estimate[2])),
    ## dispersal p = beta
    dispersal_p = rbetat(n = 10000,
                         range = quantile(demo_params_selected[["dispersal_p"]], probs = c(0.05,0.95)),
                         shape1 = param_dists$dispersal_p$Output$beta$estimate[1],
                         shape2 = param_dists$dispersal_p$Output$beta$estimate[2]),
    ## dispersal b = lognormal
    dispersal_b = rlnormt(n = 10000,
                          range = quantile(demo_params_selected[["dispersal_b"]], probs = c(0.05,0.95)),
                          meanlog = param_dists$dispersal_b$Output$lnorm$estimate[1],
                          sdlog = param_dists$dispersal_b$Output$lnorm$estimate[2]),
    ## abundance threshold = normal
    abundance_threshold = round(
      rnormt(n = 10000,
             range = quantile(demo_params_selected[["abundance_threshold"]], probs = c(0.05,0.95)),
             mean = param_dists$abundance_threshold$Output$norm$estimate[1],
             sd = param_dists$abundance_threshold$Output$norm$estimate[2])),
    ## harvest max = beta
    harvest_max = rbetat(n = 10000,
                         range = quantile(demo_params_selected[["harvest_max"]], probs = c(0.05,0.95)),
                         shape1 = param_dists$harvest_max$Output$beta$estimate[1],
                         shape2 = param_dists$harvest_max$Output$beta$estimate[2]),
    ## harvest z = norm
    harvest_z = rnormt(n = 10000,
                       range = quantile(demo_params_selected[["harvest_z"]], probs = c(0.05,0.95)),
                       mean = param_dists$harvest_z$Output$norm$estimate[1],
                       sd = param_dists$harvest_z$Output$norm$estimate[2]),
    ## humans multiplier = beta
    humans_multiplier = rbetat(n = 10000,
                               range = quantile(demo_params_selected[["humans_multiplier"]], probs = c(0.05,0.95)),
                               shape1 = param_dists$humans_multiplier$Output$beta$estimate[1],
                               shape2 = param_dists$humans_multiplier$Output$beta$estimate[2]),
    ## p = beta
    p = rbetat(n = 10000,
               range = quantile(demo_params_selected[["p"]], probs = c(0.05,0.95)),
               shape1 = param_dists$p$Output$beta$estimate[1],
               shape2 = param_dists$p$Output$beta$estimate[2])
  )}

dt_lhs_run3

summary(dt_lhs_run3)

par(mfrow = n2mfrow(ncol(dt_lhs_run3)), mar = c(2,2,2,0.5))
lapply(colnames(dt_lhs_run3), function(param) {
  orig <- demo_params_selected[[param]]
  sampled <- dt_lhs_run3[[param]]
  dOrig <- density(orig, from = min(orig), to = max(orig), adjust = 1,
                    bw = "SJ")
  dSamp <- density(sampled, from = min(sampled), to = max(sampled), adjust = 1,
                   bw = "SJ")
  p1 <- hist(orig, breaks = 20, plot = FALSE)
  p2 <- hist(sampled, breaks = 20, plot = FALSE)
  plot(p1, col = rgb(0,0,1,1/4), xlim = range(orig),
       freq = FALSE, ylim = range(p1$density, p2$density),
       main = param, xlab = NA, ylab = NA)  # first histogram
  plot(p2, col = rgb(1,0,0,1/4), xlim = range(orig), add = TRUE,
       ylim = range(p1$density, p2$density), freq= FALSE)  # second
  # plot(dOrig, ylim = c(min(c(dOrig$y, dSamp$y)), max(c(dOrig$y, dSamp$y))),
  #      main = param,
  #      col = "black", lty = 1)
  rug(orig, side = 1, col = rgb(0,0,1,1/4))
  rug(sampled, side = 3, col = rgb(1,0,0,1/4))
  lines(dOrig, col = rgb(0,0,1,1/4), lty = 2, lwd = 2)
  lines(dSamp, col = rgb(1,0,0,1/4), lty = 2, lwd = 2)
  abline(v = quantile(orig, c(0.05,0.95)), lty = 3, col = "grey70", lwd = 2)
})

fwrite(dt_lhs_run3, "../Simulations/Coelodonta/sample_params_run3_10000.csv")
