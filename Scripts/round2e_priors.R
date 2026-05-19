library(data.table)
library(here)
library(abc)
library(triangle)
library(qs2)
library(mc2d)
library(fitdistrplus)

options(scipen = 999)

abc_5params <- qs_read(here("Data/Validation/abc_round1_3params.qs2"))
round1_metrics <- fread(here("Data/Validation/round1_validation_metrics.csv"))
demo_params <- fread(here("Data/Input/sample_data_round1.csv"))

demo_params_selected <- demo_params[round1_metrics$dc, ][, `:=`(
  Region = abc_5params$region,
  Weights = 1 /
    (abc_5params$dist +
      .Machine$double.eps)
)][Region == TRUE, ]
summary(demo_params_selected)

demo_params_selected[1:10, ]

x <- demo_params_selected[["density_max"]]
w <- demo_params_selected[Region == "TRUE", ][["Weights"]] /
  sum(demo_params_selected[Region == "TRUE", ][["Weights"]])
min.max <- as.vector(quantile(x, c(0.05, 0.95)))
{
  set.seed(54)
  density_x <- density(
    x,
    from = min.max[1],
    to = min.max[2],
    bw = "SJ",
    kernel = "gaussian",
    n = 1000,
    adjust = 0.25,
    weights = w
  )
  resample_x <- approx(
    x = cumsum(density_x$y) / sum(density_x$y),
    y = density_x$x,
    n = 10000,
    yleft = min(density_x$x),
    yright = max(density_x$x),
    na.rm = TRUE
  )$y
}
summary(resample_x)

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
  F.a <- pgamma(min(range), shape = shape, rate = rate, ...)
  F.b <- pgamma(max(range), shape = shape, rate = rate, ...)
  u <- runif(n, min = F.a, max = F.b)
  qgamma(u, shape = shape, rate = rate, ...)
}

## beta
rbetat <- function(n, range, shape1, shape2, ...) {
  F.a <- pbeta(min(range), shape1, shape2, ...)
  F.b <- pbeta(max(range), shape1, shape2, ...)
  u <- runif(n, min = F.a, max = F.b)
  qbeta(u, shape1, shape2, ...)
}

## triangular
rtrianglet <- function(n, range, min, max, mode, ...) {
  F.a <- ptriang(min(range), min = min, max = max, mode = mode, ...)
  F.b <- ptriang(max(range), min = min, max = max, mode = mode, ...)
  u <- runif(n, min = F.a, max = F.b)
  qtriang(u, min = min, max = max, mode = mode, ...)
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
triD <- fitdist(
  x,
  dist = "triang",
  start = list(
    min = min(min.max),
    max = max(min.max),
    mode = raster::modal(resample_x)
  ),
  method = "mge",
  gof = "KS"
)
poisD <- fitdist(x, distr = "pois", method = "mme")
negbD <- fitdist(x, distr = "nbinom", method = "mme")

summary(lnD)
summary(triD)
summary(poisD)
summary(negbD)
cdfcomp(list(lnD, triD, negbD))
qqcomp(list(lnD, triD, negbD))
ppcomp(list(lnD, triD, negbD))
gofstat(list(lnD, triD, negbD))

## Iterate through each demo variable and check distributions
param_dists <- lapply(
  colnames(demo_params_selected)[-c(25, 31, 48, 49)],
  function(param) {
    message(param)
    x <- demo_params_selected[[param]]
    w <- demo_params_selected[Region == "TRUE", ][["Weights"]] /
      sum(demo_params_selected[Region == "TRUE", ][["Weights"]])
    min.max <- as.vector(quantile(x, c(0.05, 0.95)))
    # set seed and resample to 10,000 obs from ECDF
    ## https://stackoverflow.com/questions/32871602/r-generate-data-from-a-probability-density-distribution
    {
      set.seed(20212810)
      density_x <- density(
        x,
        from = min.max[1],
        to = min.max[2],
        bw = "SJ",
        kernel = "gaussian",
        n = 1000,
        adjust = 0.75,
        weights = w
      )
      resample_x <- approx(
        x = cumsum(density_x$y) / sum(density_x$y),
        y = density_x$x,
        n = 10000,
        yleft = min(density_x$x),
        yright = max(density_x$x),
        na.rm = TRUE
      )$y
    }
    par(mfrow = c(2, 1))
    {
      plot(density(x, from = min(x), to = max(x)), main = param)
      abline(
        v = min.max,
        col = "grey70",
        lty = 2,
        lwd = 2
      )
    }
    disc <- ifelse(
      param %in% c("density_max", "abundance_threshold"),
      TRUE,
      FALSE
    )
    if (param %in% c("density_max", "abundance_threshold")) {
      resample_x <- round(resample_x)
    }
    # Estimate parameters from interpolated/density estimated obs
    orig_dists <- list(
      normD <- tryCatch(
        fitdist(
          resample_x,
          distr = "norm",
          method = "mge",
          gof = "KS",
          discrete = disc
        ),
        error = function(err) {
          invisible({
            err
          })
        }
      ),
      lnD <- tryCatch(
        fitdist(
          resample_x,
          distr = "lnorm",
          method = "mge",
          gof = "KS",
          discrete = disc
        ),
        error = function(err) {
          invisible({
            err
          })
        }
      ),
      gD <- tryCatch(
        fitdist(
          resample_x,
          distr = "gamma",
          method = "mge",
          gof = "KS",
          discrete = disc
        ),
        error = function(err) {
          invisible({
            err
          })
        }
      ),
      bD <- tryCatch(
        fitdist(
          resample_x,
          distr = "beta",
          method = "mge",
          gof = "KS",
          discrete = disc
        ),
        error = function(err) {
          invisible({
            err
          })
        }
      ),
      triD <- tryCatch(
        fitdist(
          resample_x,
          dist = "triang",
          start = list(
            min = min(min.max),
            max = max(min.max),
            mode = raster::modal(x)
          ),
          method = "mge",
          gof = "KS",
          discrete = disc
        ),
        error = function(err) {
          invisible({
            err
          })
        }
      ),
      unifD <- tryCatch(
        fitdist(
          resample_x,
          distr = "unif",
          method = "mge",
          gof = "KS",
          discrete = disc
        ),
        error = function(err) {
          invisible({
            err
          })
        }
      ),
      poisD <- tryCatch(
        fitdist(
          resample_x,
          distr = "pois",
          method = "mme",
          discrete = TRUE
        ),
        error = function(err) {
          invisible({
            err
          })
        },
        warning = function(warn) {
          invisible({
            warn
          })
        }
      ),
      negbD <- tryCatch(
        fitdist(
          resample_x,
          distr = "nbinom",
          method = "mme",
          discrete = TRUE
        ),
        error = function(err) {
          invisible({
            err
          })
        },
        warning = function(warn) {
          invisible({
            warn
          })
        }
      )
    )
    # index for dist that errored
    keep_dists <- sapply(orig_dists, function(i) {
      inherits(i, c("error", "warning"))
    })
    orig_dists <- orig_dists[!keep_dists]
    names(orig_dists) <- c(
      "norm",
      "lnorm",
      "gamma",
      "beta",
      "triangular",
      "unif",
      "pois",
      "negbin"
    )[!keep_dists]
    cdfcomp(
      orig_dists,
      main = param,
      addlegend = TRUE,
      ylim = c(0, 1)
    )
    par(mfrow = c(1, 1))
    # gofstat(orig_dists, fitnames = names(orig_dists))
    return(
      list(Output = orig_dists)
    )
  }
)
names(param_dists) <- colnames(demo_params_selected)[-c(25, 31, 48, 49)]

# For each param, calculate the min AIC, then calc. delta AIC for each distribution
delta_aic <- lapply(names(param_dists), function(param) {
  var_dist <- param_dists[param][param][[1]]$Output
  dist_aic <- sapply(var_dist, function(i) {
    i$aic
  })
  delta_aic <- dist_aic - min(dist_aic, na.rm = TRUE)
  delta_aic <- sort(delta_aic, na.last = TRUE)
  return(delta_aic)
})
names(delta_aic) <- names(param_dists)
delta_aic

# Make 10,000 resamples based on results
{
  set.seed(2026311)

  params_to_sample <- names(delta_aic)

  integer_params <- c(
    "abundance_threshold",
    "initial_release",
    "density_max",
    "fecundity",
    "infected_t1"
  )

  sample_one_param <- function(param, n = 10000) {
    best_dist <- names(delta_aic[[param]])[1]
    est <- param_dists[[param]]$Output[[best_dist]]$estimate
    range <- as.vector(
      quantile(
        demo_params_selected[[param]],
        probs = c(0.05, 0.95),
        na.rm = TRUE
      )
    )

    draws <- switch(
      best_dist,
      norm = rnormt(
        n = n,
        range = range,
        mean = est[1],
        sd = est[2]
      ),
      lnorm = rlnormt(
        n = n,
        range = range,
        meanlog = est[1],
        sdlog = est[2]
      ),
      gamma = rgammat(
        n = n,
        range = range,
        shape = est[1],
        rate = est[2]
      ),
      beta = rbetat(
        n = n,
        range = range,
        shape1 = est[1],
        shape2 = est[2]
      ),
      triangular = rtrianglet(
        n = n,
        range = range,
        min = est[1],
        max = est[2],
        mode = est[3]
      ),
      unif = {
        lower <- max(range[1], est[1])
        upper <- min(range[2], est[2])
        if (lower >= upper) {
          lower <- range[1]
          upper <- range[2]
        }
        runif(n = n, min = lower, max = upper)
      },
      pois = rpoist(
        n = n,
        range = range,
        lambda = est[1]
      ),
      negbin = rnbt(
        n = n,
        range = range,
        size = est[1],
        mu = est[2]
      ),
      stop(paste(
        "Unsupported best-fit distribution for",
        param,
        ":",
        best_dist
      ))
    )

    if (param %chin% integer_params) {
      draws <- round(draws)
    }

    draws
  }

  sampled_list <- lapply(params_to_sample, sample_one_param, n = 10000)
  names(sampled_list) <- params_to_sample
  dt_lhs_run2 <- as.data.table(sampled_list)

  best_dist_by_param <- sapply(delta_aic, function(x) names(x)[1])
}

dt_lhs_run2

summary(dt_lhs_run2)

write.csv(
  dt_lhs_run2,
  here("Data_minimal/Input/sample_data_round2e.csv"),
  row.names = FALSE
)
