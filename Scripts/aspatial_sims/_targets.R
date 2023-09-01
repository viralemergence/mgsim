library(targets)
library(tarchetypes)
library(future)
library(here)
library(poems)
library(dplyr)
library(purrr)
library(crew)
plan(multisession)
controller <- crew::crew_controller_local(name = "my_controller",
                                          workers = 12,
                                          seconds_idle = 3)
here("Scripts/aspatial_sims/R") |> list.files(full.names = TRUE) |> sapply(source)
tar_option_set(
  garbage_collection = TRUE,
  packages = c(
    "poems",
    "dplyr",
    "fs",
    "furrr",
    "ggplot2",
    "paletteer",
    "cowplot",
    "fst",
    "purrr",
    "tibble",
    "tidyr"
  ),
  # controller = controller
)

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

sample_df <- simulation_input_df(10000, 239084)
inputs <- sample_df %>% rowwise() %>% group_split() %>% map(aspatial_siri_prep, init = init) %>% transpose()

simulations <- tar_map_rep(
  name = simulations,
  command = siri_model_year(summer_length, init, summer_params, winter_params),
  values = inputs,
  names = tidyselect::all_of("sample"),
  columns = tidyselect::all_of("sample"),
  batches = 10,
  reps = 1,
  combine = TRUE,
  format = "fst"
)

metrics <- tar_target(metrics, extract_metrics(simulations, sample_df))
validation <- tar_target(validation, validate_sims(metrics, 1e-3))
metrics_plot <- tar_target(metrics_plot, plot_validation_metrics(metrics, here(
  "Scripts/aspatial_sims/best_case_sims.RDS"
)))
selected_models_plot <- tar_target(selected_models_plot,
                                   plot_ensemble_mean(validation, folder, metrics))

list(simulations, metrics, validation, metrics_plot, selected_models_plot)
