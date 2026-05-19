library(data.table)
library(here)
library(abc)
library(triangle)
library(qs)
library(mc2d)
library(fitdistrplus)

options(scipen = 999)

demo_params <- fread(here("Data/Input/sample_data_round2b.csv"))
abc_first_pass <- qread(here("Data/Validation/abc_round2c_2.qs"))
round2_metrics <- fread(here("Data/Validation/round2c_validation_metrics.csv"))

demo_params_selected <- demo_params[round2_metrics$dc, ][, `:=`(
  Region = abc_first_pass$region,
  Weights = 1 /
    (abc_first_pass$dist +
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
  colnames(demo_params_selected)[-c(46, 47)],
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
names(param_dists) <- colnames(demo_params_selected)[-c(46, 47)]

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
  set.seed(20250116)

  dt_lhs_run2 <- data.table(
    dispersal_p_juv = rbetat(
      n = 10000,
      range = quantile(
        demo_params_selected[["dispersal_p_juv"]],
        probs = c(0.05, 0.95)
      ),
      shape1 = param_dists$dispersal_p_juv$Output$beta$estimate[1],
      shape2 = param_dists$dispersal_p_juv$Output$beta$estimate[2]
    ),
    dispersal_p_adult = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["dispersal_p_adult"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$dispersal_p_adult$Output$triangular$estimate[1],
      max = param_dists$dispersal_p_adult$Output$triangular$estimate[2],
      mode = param_dists$dispersal_p_adult$Output$triangular$estimate[3]
    ),
    dispersal_r_juv = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["dispersal_r_juv"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$dispersal_r_juv$Output$lnorm$estimate[1],
      sdlog = param_dists$dispersal_r_juv$Output$lnorm$estimate[2]
    ),
    dispersal_r_adult = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["dispersal_r_adult"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$dispersal_r_adult$Output$triangular$estimate[1],
      max = param_dists$dispersal_r_adult$Output$triangular$estimate[2],
      mode = param_dists$dispersal_r_adult$Output$triangular$estimate[3]
    ),
    dispersal_source_n_k_cutoff = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["dispersal_source_n_k_cutoff"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$dispersal_source_n_k_cutoff$Output$triangular$estimate[
        1
      ],
      max = param_dists$dispersal_source_n_k_cutoff$Output$triangular$estimate[
        2
      ],
      mode = param_dists$dispersal_source_n_k_cutoff$Output$triangular$estimate[
        3
      ]
    ),
    dispersal_source_n_k_threshold = rbetat(
      n = 10000,
      range = quantile(
        demo_params_selected[["dispersal_source_n_k_threshold"]],
        probs = c(0.05, 0.95)
      ),
      shape1 = param_dists$dispersal_source_n_k_threshold$Output$beta$estimate[
        1
      ],
      shape2 = param_dists$dispersal_source_n_k_threshold$Output$beta$estimate[
        2
      ]
    ),
    dispersal_target_n_k_cutoff = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["dispersal_target_n_k_cutoff"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$dispersal_target_n_k_cutoff$Output$triangular$estimate[
        1
      ],
      max = param_dists$dispersal_target_n_k_cutoff$Output$triangular$estimate[
        2
      ],
      mode = param_dists$dispersal_target_n_k_cutoff$Output$triangular$estimate[
        3
      ]
    ),
    dispersal_target_n_k_threshold = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["dispersal_target_n_k_threshold"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$dispersal_target_n_k_threshold$Output$triangular$estimate[
        1
      ],
      max = param_dists$dispersal_target_n_k_threshold$Output$triangular$estimate[
        2
      ],
      mode = param_dists$dispersal_target_n_k_threshold$Output$triangular$estimate[
        3
      ]
    ),
    abundance_threshold = round(
      rtrianglet(
        n = 10000,
        range = quantile(
          demo_params_selected[["abundance_threshold"]],
          probs = c(0.05, 0.95)
        ),
        min = param_dists$abundance_threshold$Output$triangular$estimate[1],
        max = param_dists$abundance_threshold$Output$triangular$estimate[2],
        mode = param_dists$abundance_threshold$Output$triangular$estimate[3]
      )
    ),
    initial_release = round(
      rtrianglet(
        n = 10000,
        range = quantile(
          demo_params_selected[["initial_release"]],
          probs = c(0.05, 0.95)
        ),
        min = param_dists$initial_release$Output$triangular$estimate[1],
        max = param_dists$initial_release$Output$triangular$estimate[2],
        mode = param_dists$initial_release$Output$triangular$estimate[3]
      )
    ),
    density_max = round(
      rtrianglet(
        n = 10000,
        range = quantile(
          demo_params_selected[["density_max"]],
          probs = c(0.05, 0.95)
        ),
        min = param_dists$density_max$Output$triangular$estimate[1],
        max = param_dists$density_max$Output$triangular$estimate[2],
        mode = param_dists$density_max$Output$triangular$estimate[3]
      )
    ),
    fecundity = round(
      rtrianglet(
        n = 10000,
        range = quantile(
          demo_params_selected[["fecundity"]],
          probs = c(0.05, 0.95)
        ),
        min = param_dists$fecundity$Output$triangular$estimate[1],
        max = param_dists$fecundity$Output$triangular$estimate[2],
        mode = param_dists$fecundity$Output$triangular$estimate[3]
      )
    ),
    infected_t1 = round(
      rtrianglet(
        n = 10000,
        range = quantile(
          demo_params_selected[["infected_t1"]],
          probs = c(0.05, 0.95)
        ),
        min = param_dists$infected_t1$Output$triangular$estimate[1],
        max = param_dists$infected_t1$Output$triangular$estimate[2],
        mode = param_dists$infected_t1$Output$triangular$estimate[3]
      )
    ),
    beta_Sa_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_Sa_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$beta_Sa_winter$Output$triangular$estimate[1],
      max = param_dists$beta_Sa_winter$Output$triangular$estimate[2],
      mode = param_dists$beta_Sa_winter$Output$triangular$estimate[3]
    ),
    beta_Sa_summer = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_Sa_summer"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$beta_Sa_summer$Output$triangular$estimate[1],
      max = param_dists$beta_Sa_summer$Output$triangular$estimate[2],
      mode = param_dists$beta_Sa_summer$Output$triangular$estimate[3]
    ),
    beta_I2_modifier = rbetat(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_I2_modifier"]],
        probs = c(0.05, 0.95)
      ),
      shape1 = param_dists$beta_I2_modifier$Output$beta$estimate[1],
      shape2 = param_dists$beta_I2_modifier$Output$beta$estimate[2]
    ),
    mortality_Sj_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_Sj_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_Sj_winter$Output$triangular$estimate[1],
      max = param_dists$mortality_Sj_winter$Output$triangular$estimate[2],
      mode = param_dists$mortality_Sj_winter$Output$triangular$estimate[3]
    ),
    mortality_Sa_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_Sa_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_Sa_winter$Output$triangular$estimate[1],
      max = param_dists$mortality_Sa_winter$Output$triangular$estimate[2],
      mode = param_dists$mortality_Sa_winter$Output$triangular$estimate[3]
    ),
    mortality_Sj_summer = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_Sj_summer"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_Sj_summer$Output$triangular$estimate[1],
      max = param_dists$mortality_Sj_summer$Output$triangular$estimate[2],
      mode = param_dists$mortality_Sj_summer$Output$triangular$estimate[3]
    ),
    mortality_I1j_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I1j_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_I1j_winter$Output$triangular$estimate[1],
      max = param_dists$mortality_I1j_winter$Output$triangular$estimate[2],
      mode = param_dists$mortality_I1j_winter$Output$triangular$estimate[3]
    ),
    mortality_I1j_summer = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I1j_summer"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$mortality_I1j_summer$Output$lnorm$estimate[1],
      sdlog = param_dists$mortality_I1j_summer$Output$lnorm$estimate[2]
    ),
    mortality_I1a_summer = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I1a_summer"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$mortality_I1a_summer$Output$lnorm$estimate[1],
      sdlog = param_dists$mortality_I1a_summer$Output$lnorm$estimate[2]
    ),
    mortality_I1a_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I1a_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_I1a_winter$Output$triangular$estimate[1],
      max = param_dists$mortality_I1a_winter$Output$triangular$estimate[2],
      mode = param_dists$mortality_I1a_winter$Output$triangular$estimate[3]
    ),
    mortality_I2_modifier = rbetat(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I2_modifier"]],
        probs = c(0.05, 0.95)
      ),
      shape1 = param_dists$mortality_I2_modifier$Output$beta$estimate[1],
      shape2 = param_dists$mortality_I2_modifier$Output$beta$estimate[2]
    ),
    mortality_I2j_summer = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I2j_summer"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$mortality_I2j_summer$Output$lnorm$estimate[1],
      sdlog = param_dists$mortality_I2j_summer$Output$lnorm$estimate[2]
    ),
    mortality_I2j_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I2j_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_I2j_winter$Output$triangular$estimate[1],
      max = param_dists$mortality_I2j_winter$Output$triangular$estimate[2],
      mode = param_dists$mortality_I2j_winter$Output$triangular$estimate[3]
    ),
    mortality_I2a_winter = rgammat(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I2a_winter"]],
        probs = c(0.05, 0.95)
      ),
      shape = param_dists$mortality_I2a_winter$Output$gamma$estimate[1],
      rate = param_dists$mortality_I2a_winter$Output$gamma$estimate[2]
    ),
    mortality_I2a_summer = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_I2a_summer"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$mortality_I2a_summer$Output$lnorm$estimate[1],
      sdlog = param_dists$mortality_I2a_summer$Output$lnorm$estimate[2]
    ),
    mortality_Rj_summer = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_Rj_summer"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_Rj_summer$Output$triangular$estimate[1],
      max = param_dists$mortality_Rj_summer$Output$triangular$estimate[2],
      mode = param_dists$mortality_Rj_summer$Output$triangular$estimate[3]
    ),
    mortality_Rj_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_Rj_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_Rj_winter$Output$triangular$estimate[1],
      max = param_dists$mortality_Rj_winter$Output$triangular$estimate[2],
      mode = param_dists$mortality_Rj_winter$Output$triangular$estimate[3]
    ),
    mortality_Ra_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["mortality_Ra_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$mortality_Ra_winter$Output$triangular$estimate[1],
      max = param_dists$mortality_Ra_winter$Output$triangular$estimate[2],
      mode = param_dists$mortality_Ra_winter$Output$triangular$estimate[3]
    ),
    beta_Sj_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_Rj_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$beta_Rj_winter$Output$triangular$estimate[1],
      max = param_dists$beta_Rj_winter$Output$triangular$estimate[2],
      mode = param_dists$beta_Rj_winter$Output$triangular$estimate[3]
    ),
    beta_Sj_summer = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_Sj_summer"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$beta_Sj_summer$Output$triangular$estimate[1],
      max = param_dists$beta_Sj_summer$Output$triangular$estimate[2],
      mode = param_dists$beta_Sj_summer$Output$triangular$estimate[3]
    ),
    beta_Ra_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_Ra_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$beta_Ra_winter$Output$triangular$estimate[1],
      max = param_dists$beta_Ra_winter$Output$triangular$estimate[2],
      mode = param_dists$beta_Ra_winter$Output$triangular$estimate[3]
    ),
    beta_Rj_winter = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_Rj_winter"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$beta_Rj_winter$Output$lnorm$estimate[1],
      sdlog = param_dists$beta_Rj_winter$Output$lnorm$estimate[2]
    ),
    beta_Ra_summer = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_Ra_summer"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$beta_Ra_summer$Output$lnorm$estimate[1],
      sdlog = param_dists$beta_Ra_summer$Output$lnorm$estimate[2]
    ),
    beta_Rj_summer = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["beta_Rj_summer"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$beta_Rj_summer$Output$lnorm$estimate[1],
      sdlog = param_dists$beta_Rj_summer$Output$lnorm$estimate[2]
    ),
    recovery_I1j_summer = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["recovery_I1j_summer"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$recovery_I1j_summer$Output$triangular$estimate[1],
      max = param_dists$recovery_I1j_summer$Output$triangular$estimate[2],
      mode = param_dists$recovery_I1j_summer$Output$triangular$estimate[3]
    ),
    recovery_I1a_summer = rlnormt(
      n = 10000,
      range = quantile(
        demo_params_selected[["recovery_I1a_summer"]],
        probs = c(0.05, 0.95)
      ),
      meanlog = param_dists$recovery_I1a_summer$Output$lnorm$estimate[1],
      sdlog = param_dists$recovery_I1a_summer$Output$lnorm$estimate[2]
    ),
    recovery_I2j_summer = rbetat(
      n = 10000,
      range = quantile(
        demo_params_selected[["recovery_I2j_summer"]],
        probs = c(0.05, 0.95)
      ),
      shape1 = param_dists$recovery_I2j_summer$Output$beta$estimate[1],
      shape2 = param_dists$recovery_I2j_summer$Output$beta$estimate[2]
    ),
    recovery_I2a_summer = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["recovery_I1j_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$recovery_I1j_winter$Output$triangular$estimate[1],
      max = param_dists$recovery_I1j_winter$Output$triangular$estimate[2],
      mode = param_dists$recovery_I1j_winter$Output$triangular$estimate[3]
    ),
    recovery_I1j_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["recovery_I1j_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$recovery_I1j_winter$Output$triangular$estimate[1],
      max = param_dists$recovery_I1j_winter$Output$triangular$estimate[2],
      mode = param_dists$recovery_I1j_winter$Output$triangular$estimate[3]
    ),
    recovery_I1a_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["recovery_I1a_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$recovery_I1a_winter$Output$triangular$estimate[1],
      max = param_dists$recovery_I1a_winter$Output$triangular$estimate[2],
      mode = param_dists$recovery_I1a_winter$Output$triangular$estimate[3]
    ),
    recovery_I2j_winter = rtrianglet(
      n = 10000,
      range = quantile(
        demo_params_selected[["recovery_I2j_winter"]],
        probs = c(0.05, 0.95)
      ),
      min = param_dists$recovery_I2j_winter$Output$triangular$estimate[1],
      max = param_dists$recovery_I2j_winter$Output$triangular$estimate[2],
      mode = param_dists$recovery_I2j_winter$Output$triangular$estimate[3]
    ),
    recovery_I2a_winter = rgammat(
      n = 10000,
      range = quantile(
        demo_params_selected[["recovery_I2a_winter"]],
        probs = c(0.05, 0.95)
      ),
      shape = param_dists$recovery_I2a_winter$Output$gamma$estimate[1],
      rate = param_dists$recovery_I2a_winter$Output$gamma$estimate[2]
    )
  )
}

dt_lhs_run2

summary(dt_lhs_run2)

write.csv(
  dt_lhs_run2,
  here("Data_minimal/Input/sample_data_round3a.csv"),
  row.names = FALSE
)
