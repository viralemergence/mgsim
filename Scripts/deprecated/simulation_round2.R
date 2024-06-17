library(raster)
library(terra)
library(tidyverse)
library(poems)
library(epizootic)
library(qs)
library(here)
data_dir <- here("Data/Input")
parallel_cores <- 18
nsims <- 3211
burn_in_steps <- 0
timesteps <- 23 + burn_in_steps
random_seed <- 324
results_dir <- here("Data/Output/Round 2")
region <- data_dir %>% file.path("finch_region.qs") %>% qread()
env_corr <- SpatialCorrelation$new(region = region,
                                   amplitude = 0.99,
                                   breadth = 850,
                                   distance_scale = 1000)
env_corr$calculate_compact_decomposition(decimals = 4)
# burn_in <- stack(replicate(burn_in_steps, raster(
#     file.path(data_dir, "breeding_season_length.tif"), layer = 1
#   )))
bsl_raster <- data_dir %>%
  file.path("breeding_season_length.tif") %>%
  stack(bands = 55:77) %>%
  #  stack(burn_in, .) |>
  mask(region$region_raster) |>
  calc(fun = round)
model_template <- DiseaseModel$new(
  simulation_function = "disease_simulator",
  region = region,
  time_steps = timesteps,
  populations = region$region_cells,
  replicates = 1,
  stages = 2,
  compartments = 4,
  seasons = 2,
  mortality_unit = list(c(1, 1, 0, 0, 1, 1, 0, 0),
                        c(1, 1, 0, 0, 1, 1, 0, 0)),
  fecundity_unit = rep(1, 8),
  fecundity_mask = rep(c(0, 1), 4),
  transmission_unit = rep(0, 4),
  transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0),
  recovery_unit = rep(0, 4),
  recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1),
  breeding_season_length = bsl_raster,
  simulation_order = list(c("dispersal", "season_functions", "results"),
                          c("transition", "season_functions", "results")),
  results_selection = c("abundance"),
  results_breakdown = "segments",
  season_functions = list(siri_model_winter, siri_model_summer),
  hs_file = "habitat_suitability_v3",
  mask_file = "native_range_mask",
  verbose = FALSE,
  attribute_aliases = list(
    mortality_Sj_summer = "mortality$summer$a",
    mortality_Sa_summer = "mortality$summer$b",
    mortality_I1j_summer = "mortality$summer$c",
    mortality_I1a_summer = "mortality$summer$d",
    mortality_Rj_summer = "mortality$summer$e",
    mortality_Ra_summer = "mortality$summer$f",
    mortality_I2j_summer = "mortality$summer$g",
    mortality_I2a_summer = "mortality$summer$h",
    mortality_Sj_winter = "mortality$nonbreeding$a",
    mortality_Sa_winter = "mortality$nonbreeding$b",
    mortality_I1j_winter = "mortality$nonbreeding$c",
    mortality_I1a_winter = "mortality$nonbreeding$d",
    mortality_Rj_winter = "mortality$nonbreeding$e",
    mortality_Ra_winter = "mortality$nonbreeding$f",
    mortality_I2j_winter = "mortality$nonbreeding$g",
    mortality_I2a_winter = "mortality$nonbreeding$h",
    beta_Sj_summer = "transmission$summer$a",
    beta_Sa_summer = "transmission$summer$b",
    beta_Rj_summer = "transmission$summer$c",
    beta_Ra_summer = "transmission$summer$d",
    beta_Sj_winter = "transmission$nonbreeding$a",
    beta_Sa_winter = "transmission$nonbreeding$b",
    beta_Rj_winter = "transmission$nonbreeding$c",
    beta_Ra_winter = "transmission$nonbreeding$d",
    recovery_I1j_summer = "recovery$summer$a",
    recovery_I1a_summer = "recovery$summer$b",
    recovery_I2j_summer = "recovery$summer$c",
    recovery_I2a_summer = "recovery$summer$d",
    recovery_I1j_winter = "recovery$nonbreeding$a",
    recovery_I1a_winter = "recovery$nonbreeding$b",
    recovery_I2j_winter = "recovery$nonbreeding$c",
    recovery_I2a_winter = "recovery$nonbreeding$d",
    dispersal1 = "dispersal$a",
    dispersal2 = "dispersal$b"
  )
)
b_lookup <- data.frame(d_max = -Inf, b = 0:904)
for (i in 2:904) {
  b_lookup$d_max[i] <- which.max(exp(-1*(1:1501)/b_lookup$b[i]) <= 0.19)
}
b_lookup$d_max[905] <- 1501

# distance_matrix <- dispersal_gen$calculate_distance_matrix()
# dispersal_gen$calculate_distance_data(distance_matrix = distance_matrix)
# these are pre-calculated because they're computationally intensive
distance_data <- qread(file.path(data_dir, "dispersal_distance_data.qs"))

adult_dispersal_gen <- DispersalGenerator$new(
  region = region,
  spatial_correlation = env_corr,
  distance_classes = seq(10, 1500, 10),
  distance_scale = 1000, # km
  dispersal_function_data = b_lookup,
  inputs = c("dispersal_p_adult",
             "dispersal_r_adult"),
  attribute_aliases = list(dispersal_r_adult = "dispersal_max_distance",
                           dispersal_p_adult = "dispersal_proportion",
                           dispersal_adult = "dispersal_data"),
  decimals = 3
)
adult_dispersal_gen$distance_data <- distance_data

juvenile_dispersal_gen <- DispersalGenerator$new(
  region = region,
  spatial_correlation = env_corr,
  distance_classes = seq(10, 1500, 10),
  distance_scale = 1000, # km
  dispersal_function_data = b_lookup,
  decimals = 3,
  inputs = c("dispersal_p_juv",
             "dispersal_r_juv"),
  attribute_aliases = list(dispersal_r_juv = "dispersal_max_distance",
                           dispersal_p_juv = "dispersal_proportion",
                           dispersal_source_n_k_cutoff = "dispersal_source_n_k$cutoff",
                           dispersal_juv = "dispersal_data"),
  decimals = 3
)
juvenile_dispersal_gen$distance_data <- distance_data

# Test run
adult_dispersal_data <- adult_dispersal_gen$generate(
  input_values = list(
    dispersal_p_adult = 0.5,
    dispersal_r_adult = 300
  )
) |> _$dispersal_data
head(adult_dispersal_data[[1]])
capacity_gen <- Generator$new(
  description = "capacity",
  region = region,
  generate_rasters = FALSE,
  # use but don't generate
  burn_in_steps = burn_in_steps,
  round1_dir = here::here("Data/Output/Round 1"),
  time_steps = timesteps,
  generative_requirements = list(
    hs_raster = "file",
    round1_abundance = "file",
    initial_abundance = "function",
    carrying_capacity = "function"
  ),
  inputs = c("density_max", "hs_file", "abundance_file", "infected_t1"),
  outputs = c("initial_abundance", "carrying_capacity")
)
# Here we tell the generator to import the HS file and save it as "hs_matrix"
capacity_gen$add_file_template(
  param = "hs_raster",
  path_template = file.path(data_dir, "%s.tif"),
  path_params = "hs_file",
  file_type = "tif"
)
capacity_gen$add_file_template(
  param = "round1_abundance",
  path_template = file.path(capacity_gen$get_attribute("round1_dir"),
                            "sample_%s_results.qs"),
  path_params = "abundance_file",
  file_type = "QS"
)
# Here we subset the hs_matrix to have only the region cells, and we add the burn in
# Also, we tell the generator to generate the carrying_capacity based on "density_max" and "hs_matrix".
capacity_gen$add_function_template(
  param = "carrying_capacity",
  function_def = function(params) {
    hs_matrix <- params$hs_raster |> as.matrix() |>
      _[params$region$region_indices, 55:(54 + params$time_steps - params$burn_in_steps)]
    hs_matrix[!is.finite(hs_matrix)] <- 0
    # repeat the first timestep n times as burn in
    if (params$burn_in_steps > 1) {
      hs_matrix <- cbind(replicate(params$burn_in_steps, hs_matrix[, 1]), hs_matrix)
    }
    # round the density values
    round(params$density_max * hs_matrix)
  },
  call_params = c("density_max", "hs_raster", "burn_in_steps", "region",
                  "time_steps")
)
# Here we tell the generator what function to use to generate initial_abundance
# based on the carrying capacity of the first time step, the native range, and
# the number of finches released in Jones Beach, NY.
capacity_gen$add_function_template(
  param = "initial_abundance",
  function_def = function(params) {
    stages <- params$round1_abundance$abundance_stages
    adults <- stages$stage_2[ , ncol(stages$stage_2), 2]
    juv <- stages$stage_1[ , ncol(stages$stage_1), 2]
    infected_adults <- round(runif(1, 0, params$infected_t1))
    infected_juv <- params$infected_t1 - infected_adults
    initial_matrix <- matrix(0, nrow = 8, ncol = 6355)
    initial_matrix[2, ] <- adults
    initial_matrix[1, ] <- juv
    initial_matrix[3, 3531] <- infected_juv
    initial_matrix[4, 3531] <- infected_adults
    return(initial_matrix)
  },
  call_params = c("round1_abundance", "infected_t1")
)

system.time({
  test_capacity <- capacity_gen$generate(input_values = list(density_max = 186000,
                                                             infected_t1 = 5,
                                                             hs_file = "habitat_suitability_v3",
                                                             abundance_file = 4615))
})

raster::plot(
  region$raster_from_values(test_capacity[[1]][1,]),
  main = "Initial abundance of susceptible juveniles",
  colNA = "blue"
)
raster::plot(
  region$raster_from_values(test_capacity[[1]][3,]),
  main = "Initial abundance of infected juveniles",
  colNA = "blue"
)
lhs_generator <- LatinHypercubeSampler$new()

# Transmission parameters
lhs_generator$set_uniform_parameter("beta_Sa_winter", lower = 0, upper = 0.07588)
lhs_generator$set_uniform_parameter("beta_Sa_summer", lower = 0, upper = 0.007784)
lhs_generator$set_triangular_parameter("Sj_multiplier", lower = 0, upper = 8.5,
                                       mode = 3)
lhs_generator$set_beta_parameter("beta_I2_modifier", alpha = 1.547023,
                                 beta = 0.4239236)

# Recovery parameters
lhs_generator$set_beta_parameter("recovery_I1", alpha = 9.347533,
                                 beta = 620.1732)
lhs_generator$set_beta_parameter("recovery_I2", alpha = 1.181112,
                                 beta = 29.18489)

# How many birds are infected in DC at timestep 1?
lhs_generator$set_uniform_parameter("infected_t1", lower = 1, upper = 20, decimals = 0)

round1_sample_data <- read_csv(file.path(data_dir, "sample_data_round1.csv")) |>
  left_join(read_csv(here("Data/Validation/dc_round1.csv")),
            by = c("sample" = "sim")) |>
  filter(dc)

sample_data <- lhs_generator$generate_samples(number = nsims,
                                              random_seed = random_seed) |>
  mutate(beta_Sj_winter = beta_Sa_winter * Sj_multiplier,
         beta_Sj_summer = beta_Sa_summer * Sj_multiplier,
         beta_Ra_winter = beta_Sa_winter * beta_I2_modifier,
         beta_Rj_winter = beta_Sj_winter * beta_I2_modifier,
         beta_Ra_summer = beta_Sa_summer * beta_I2_modifier,
         beta_Rj_summer = beta_Sj_summer * beta_I2_modifier,
         recovery_I1j_summer = recovery_I1,
         recovery_I1a_summer = recovery_I1,
         recovery_I2j_summer = recovery_I2,
         recovery_I2a_summer = recovery_I2,
         recovery_I1j_winter = recovery_I1,
         recovery_I1a_winter = recovery_I1,
         recovery_I2j_winter = recovery_I2,
         recovery_I2a_winter = recovery_I2) |>
  bind_cols(round1_sample_data) |>
  rename(abundance_file = sample)
write_csv(sample_data, file.path(data_dir, "sample_data_round2.csv"))
handler <- SimulationHandler$new(
  sample_data = sample_data,
  model_template = model_template,
  generators = list(juvenile_dispersal_gen,
                    adult_dispersal_gen,
                    capacity_gen),
  parallel_cores = parallel_cores,
  results_dir = results_dir
)
sim_log <- handler$run()
