library(tidyverse)
library(qs)
library(ks)

demo_params <- read_csv(here("Data/Input/sample_data_round2b.csv"))
abc_first_pass <- qread(here("Data/Validation/abc_round2c_2.qs"))
round2_metrics <- read_csv(here(
  "Data/Validation/round2c_validation_metrics.csv"
))

demo_params_selected <- demo_params[round2_metrics$dc, ] |>
  mutate(
    Region = abc_first_pass$region,
    Weights = 1 /
      (abc_first_pass$dist +
        .Machine$double.eps)
  ) |>
  filter(Region == TRUE)

# KDE-based sampling approach

set.seed(123) # for reproducibility
n_samples <- 10000

# Function to sample from KDE for a single parameter
sample_kde <- function(x, weights, n) {
  # Fit weighted KDE using ks package
  normalized_weights <- weights / sum(weights)

  # Calculate weighted mean and standard deviation for bandwidth
  weighted_mean <- sum(x * normalized_weights)
  weighted_var <- sum(normalized_weights * (x - weighted_mean)^2)
  weighted_sd <- sqrt(weighted_var)

  # Use Silverman's rule of thumb for bandwidth
  h <- 1.06 * weighted_sd * length(x)^(-1 / 5)

  # Create fine grid for KDE evaluation
  grid <- seq(min(x) - 3 * h, max(x) + 3 * h, length.out = 512)

  # Evaluate weighted KDE
  kde_vals <- sapply(grid, function(xi) {
    sum(normalized_weights * dnorm((x - xi) / h)) / h
  })

  # Normalize to ensure it's a proper probability distribution
  kde_vals <- kde_vals / sum(kde_vals)

  # Sample from the fitted distribution
  sample(grid, size = n, replace = TRUE, prob = kde_vals)
}

round3b_priors <- demo_params_selected |>
  select(-Region, -Weights) |>
  reframe(across(
    everything(),
    ~ sample_kde(., demo_params_selected$Weights, n_samples)
  ))

round3b_priors |>
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Value"
  ) |>
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~Parameter, scales = "free") +
  theme_minimal()

round3b_priors |>
  write_csv(here("Data/Input/sample_data_round3b.csv"))
