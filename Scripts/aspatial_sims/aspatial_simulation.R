library(here)
library(poems)
library(furrr)
library(tidyverse)
library(paletteer)
library(cowplot)
library(fst)
library(fs)
plan(multisession)
here("Scripts/aspatial_sims/R") |> list.files(full.names = TRUE) |> sapply(source)
folder <- "/Users/caryinstitute/Documents/Very_Large_Data/aspatial_sims"
init <- c(
  Sa = 2000,
  Sj = 0,
  I1j = 0,
  I1a = 1,
  Rj = 0,
  Ra = 0,
  I2j = 0,
  I2a = 0
)
sample_data <- simulation_input_df(10000, 98234)
inputs <- sample_data %>% rowwise() %>% group_split() %>% 
  map(aspatial_siri_prep, init = init)
aspatial_sim(inputs, 10000, 50, folder)
metrics <- extract_metrics(folder, sample_data, "cases")
validation <- validate_sims(metrics, 2e-4)
metrics_plot <- plot_validation_metrics(metrics, here("Scripts/aspatial_sims/best_case_sims.RDS"))
ensemble_plot <- plot_ensemble_mean(validation, folder, metrics)
top_models_plot <- plot_top_models(validation, folder, metrics, models = 25)