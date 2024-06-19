#### Setup ####
library(tidyverse)
library(metaRange)
library(epizootic)
library(poems)
library(qs)
library(terra)
library(tictoc)
data_dir <- "/Users/caryinstitute/Documents/mgsim/Data/Input"
results_dir <- here::here("Data/Output/epizootic_test")
random_seed <- 90
n_sims <- 10000
region <- data_dir |> file.path("finch_region.qs") |> qread()

#### Create sample data frame ####
lhs_generator <- LatinHypercubeSampler$new()

# Dispersal parameters
lhs_generator$set_beta_parameter("dispersal_p_juv", alpha = 9.834837, 
                                 beta = 2.019125)
lhs_generator$set_beta_parameter("dispersal_p_adult", alpha = 1.5685, 
                                 beta = 2.365266)
lhs_generator$set_truncnorm_parameter("dispersal_r_juv", lower = 0, upper = 1500, 
                                      mean = 725.9071, sd = sqrt(204006.6))
lhs_generator$set_normal_parameter("dispersal_r_adult", mean = 679.4172,
                                   sd = sqrt(18594.59))
lhs_generator$set_uniform_parameter("dispersal_source_n_k_cutoff", lower = 0, 
                                    upper = 1)
lhs_generator$set_uniform_parameter("dispersal_source_n_k_threshold", lower = 0, 
                                    upper = 1)
lhs_generator$set_uniform_parameter("dispersal_target_n_k_cutoff", lower = 0, 
                                    upper = 1)
lhs_generator$set_uniform_parameter("dispersal_target_n_k_threshold", lower = 0, 
                                    upper = 1)

# Population growth parameters
lhs_generator$set_uniform_parameter("abundance_threshold", lower = 0, upper = 25, decimals = 0)
lhs_generator$set_uniform_parameter("initial_release", lower = 5, upper = 50, decimals = 0)
lhs_generator$set_uniform_parameter("density_max", lower = 186000, upper = 310000, decimals = 0)
lhs_generator$set_poisson_parameter("fecundity", lambda = 8.509018)

# How many birds are infected in DC at timestep 1?
lhs_generator$set_uniform_parameter("infected_t1", lower = 1, upper = 20, decimals = 0)

# Transmission parameters
lhs_generator$set_uniform_parameter("beta_Sa_winter", lower = 0, upper = 0.07588)
lhs_generator$set_uniform_parameter("beta_Sa_summer", lower = 0, upper = 0.007784)
lhs_generator$set_triangular_parameter("Sj_multiplier", lower = 0, upper = 8.5,
                                       mode = 3)
lhs_generator$set_beta_parameter("beta_I2_modifier", alpha = 1.547023,
                                 beta = 0.4239236)

# Mortality parameters
lhs_generator$set_beta_parameter("mortality_Sj_winter", alpha = 3.962104,
                                 beta = 2.228683)
lhs_generator$set_beta_parameter("mortality_Sa_winter", alpha = 21.89136,
                                 beta = 19.59278)
lhs_generator$set_beta_parameter("mortality_Sj_summer", alpha = 14.51403,
                                 beta = 21.53632)
lhs_generator$set_beta_parameter("mortality_I1j_summer", alpha = 2.756404,
                                 beta = 62.47181)
lhs_generator$set_beta_parameter("mortality_I1j_winter", alpha = 2.756404,
                                 beta = 62.47181)
lhs_generator$set_beta_parameter("mortality_I1a_summer", alpha = 1.771183,
                                 beta = 27.19457)
lhs_generator$set_beta_parameter("mortality_I1a_winter", alpha = 1.678424,
                                 beta = 41.15975)
lhs_generator$set_beta_parameter("mortality_I2_modifier", alpha = 1.033367,
                                 beta = 3.505319)

# Recovery parameters
lhs_generator$set_beta_parameter("recovery_I1", alpha = 9.347533,
                                 beta = 620.1732)
lhs_generator$set_beta_parameter("recovery_I2", alpha = 1.181112,
                                 beta = 29.18489)

sample_data <- lhs_generator$generate_samples(number = n_sims, 
                                              random_seed = random_seed) |>
  mutate(mortality_Sa_summer = 0, 
        mortality_I2j_summer = mortality_I2_modifier * mortality_I1j_summer,
        mortality_I2j_winter = mortality_I2_modifier * mortality_I1j_winter,
        mortality_I2a_winter = mortality_I2_modifier * mortality_I1a_winter,
        mortality_I2a_summer = mortality_I2_modifier * mortality_I1a_summer,
        mortality_Rj_summer = mortality_Sj_summer,
        mortality_Ra_summer = mortality_Sa_summer,
        mortality_Rj_winter = mortality_Sj_winter,
        mortality_Ra_winter = mortality_Sa_winter,
        beta_Sj_winter = beta_Sa_winter * Sj_multiplier,
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
    select(-c("Sj_multiplier", "recovery_I1", "recovery_I2"))

write_csv(sample_data, file.path(data_dir, "sample_data_round1.csv"))

#### Create initial abundance generator ####
abundance_gen <- Generator$new(
  description = "abundance",
  region = region,
  generate_rasters = FALSE,
  burn_in_steps = 5,
  time_steps = 82,
  generative_requirements = list(
    hs_raster = "file",
    mask_raster = "file",
    Sa_abundance = "function"
  ),
  hs_file = "habitat_suitability_v3",
  mask_file = "native_range_mask",
  inputs = c("density_max", "initial_release"),
  outputs = c("Sa_abundance")
)
# Here we tell the generator to import the HS file and save it as "hs_matrix"
abundance_gen$add_file_template(
  param = "hs_raster",
  path_template = file.path(data_dir, "%s.tif"),
  path_params = "hs_file",
  file_type = "tif"
)
abundance_gen$add_file_template(
  param = "mask_raster",
  path_template = file.path(data_dir, "%s.tif"),
  path_params = "mask_file",
  file_type = "tif"
)
# Here we tell the generator what function to use to generate initial_abundance
# based on the carrying capacity of the first time step, the native range, and 
# the number of finches released in Jones Beach, NY.
abundance_gen$add_function_template(
  param = "Sa_abundance",
  function_def = function(params) {
    hs_matrix <- params$hs_raster |> 
                 raster::mask(params$mask_raster, updatevalue = 0) |> 
                 as.array() |>
                 _[,,1] |>
                 base::`*`(params$density_max) |>
                 round()
    hs_matrix[38, 118] <- params$initial_release
    return(hs_matrix)
  },
  call_params = c("hs_raster", "mask_raster", "density_max", "initial_release")
)

system.time({
  test_capacity <- abundance_gen$generate(input_values = list(density_max = 186000, 
                                                             initial_release = 50))
})

plot(
  rast(test_capacity[[1]])
)

#### Create dispersal generators ####
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
  distance_classes = seq(10, 1500, 10),
  distance_scale = 1000, # km
  dispersal_function_data = b_lookup,
  decimals = 3,
  inputs = c("dispersal_p_adult",
             "dispersal_r_adult"),
  attribute_aliases = list(dispersal_r_adult = "dispersal_max_distance",
                           dispersal_p_adult = "dispersal_proportion",
                           dispersal_adult = "dispersal_data")
)
adult_dispersal_gen$distance_data <- distance_data

juvenile_dispersal_gen <- DispersalGenerator$new(
  region = region,
  distance_classes = seq(10, 1500, 10),
  distance_scale = 1000, # km
  dispersal_function_data = b_lookup,
  decimals = 3,
  inputs = c("dispersal_p_juv",
             "dispersal_r_juv"),
  attribute_aliases = list(dispersal_r_juv = "dispersal_max_distance",
                 dispersal_p_juv = "dispersal_proportion",
                 dispersal_source_n_k_cutoff = "dispersal_source_n_k$cutoff",
                 dispersal_juv = "dispersal_data")
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

#### Create simulation object ####
# Creates the environment 
landscape <- c(file.path(data_dir, "breeding_season_length.tif"), 
  file.path(data_dir, "habitat_suitability_v3.tif")) |> map(rast) |> 
  map(\(x) x[[1:77]]) |> map_at(1, round) |> sds()
names(landscape) <- c("breeding_season_length", "habitat_suitability")
# Create the object
sim <- create_simulation(
  source_environment = landscape,
  ID = "test_simulation",
  seed = 1566
)
sim$set_time_layer_mapping(c(rep_len(1, 5), 2:77))
sim$add_species("house_finch")

#### Add invariant traits ####
# Here I add traits that will be the same across all simulations
r <- region$region_raster
row_indices <- rowFromCell(r, region$region_indices)
col_indices <- colFromCell(r, region$region_indices)
index_matrix <- cbind(row_indices, col_indices)

sim$add_traits(
  species = "house_finch",
  population_level = TRUE,
  Sj_abundance = 0,
  # Sa_abundance is missing because it varies across simulations
  I1j_abundance = 0,
  I1a_abundance = 0,
  Rj_abundance = 0,
  Ra_abundance = 0,
  I2j_abundance = 0,
  I2a_abundance = 0
)

sim$add_traits(
  species = "house_finch",
  population_level = FALSE,
  index_matrix = index_matrix
)

#### Add one-time processes ####
# Here I add processes for disease introduction, which happens only once
sim$add_process(
  species = "house_finch", 
  process_name = "disease_introduction",
  process_fun = function() {
    self$traits$Sa_abundance[44, 113] <- self$traits$Sa_abundance[44, 113] - self$traits$infected_t1
    self$traits$I1a_abundance[44, 113] <- self$traits$infected_t1
  },
  execution_priority = 2,
  queue = FALSE
)
sim$add_process(
  process_name = "activate_disease",
  process_fun = function() {
    if (self$get_current_time_step() == 54) {
      if (self$house_finch$traits$Sa_abundance[44, 113] >= self$house_finch$traits$infected_t1) {
        self$queue$enqueue(self$house_finch$processes$disease_introduction)
      } else {
        self$exit()
      }
    }
  },
  execution_priority = 1
)
sim$add_process(
  process_name = "stop_disease_introduction",
  process_fun = function() {
    if (self$get_current_time_step() == 56) {
      self$queue$dequeue(self$house_finch$processes$disease_introduction$get_PID())
    }
  },
  execution_priority = 1
)

#### Add results-saving processes ####
sim$add_process(
  process_name = "save_abundance_winter",
  process_fun = function() {
    if (!dir.exists(self$globals$results_dir)) {
      dir.create(self$globals$results_dir, recursive = TRUE)
    }

    if (!dir.exists(self$globals$results_dir)) {
      stop("The results directory could not be created.")
    }

    if (self$get_current_time_step() > 54) {
      save_species(
        x = self$house_finch,
        traits = c("Sj_abundance", "Sa_abundance", "I1j_abundance", "I1a_abundance", 
                  "Rj_abundance", "Ra_abundance", "I2j_abundance", "I2a_abundance"),
        prefix = paste0("winter_", self$get_current_time_step(), "_"),
        path = self$globals$results_dir,
        overwrite = TRUE
      )
    } else {
      save_species(
        x = self$house_finch,
        traits = c("Sj_abundance", "Sa_abundance"),
        prefix = paste0("winter_", self$get_current_time_step(), "_"),
        path = self$globals$results_dir,
        overwrite = TRUE
      )
    }
  },
  execution_priority = 5
)
sim$add_process(
  process_name = "save_abundance_summer",
  process_fun = function() {

    if (self$get_current_time_step() > 54) {
      save_species(
        x = self$house_finch,
        traits = c("Sj_abundance", "Sa_abundance", "I1j_abundance", "I1a_abundance", 
                  "Rj_abundance", "Ra_abundance", "I2j_abundance", "I2a_abundance"),
        prefix = paste0("summer_", self$get_current_time_step(), "_"),
        path = self$globals$results_dir,
        overwrite = TRUE
      )
    } else {
      save_species(
        x = self$house_finch,
        traits = c("Sj_abundance", "Sa_abundance"),
        prefix = paste0("summer_", self$get_current_time_step(), "_"),
        path = self$globals$results_dir,
        overwrite = TRUE
      )
    }
  },
  execution_priority = 8
)

#### Add dispersal process ####
sim$add_process(
  species = "house_finch", 
  process_name = "dispersal",
  process_fun = function() {
    dispersal_fun <- disease_dispersal(
      replicates = 1,
      time_steps = 77,
      populations = 6355,
      demographic_stochasticity = TRUE,
      dispersal = list(self$traits$dispersal1,
                       self$traits$dispersal2),
      dispersal_type = "stages",
      dispersal_source_n_k = list(cutoff = self$traits$dispersal_source_n_k_cutoff,
                                  threshold = self$traits$dispersal_source_n_k_threshold),
      dispersal_target_k = NULL,
      dispersal_target_n = NULL,
      dispersal_target_n_k = list(cutoff = self$traits$dispersal_target_n_k_cutoff,
                                  threshold = self$traits$dispersal_target_n_k_threshold),
      stages = 2,
      compartments = 4,
      simulator = SimulatorReference$new()
    )

    segment_abundance <- matrix(c(self$traits$Sj_abundance[self$traits$index_matrix],
                    self$traits$Sa_abundance[self$traits$index_matrix],
                    self$traits$I1j_abundance[self$traits$index_matrix],
                    self$traits$I1a_abundance[self$traits$index_matrix],
                    self$traits$Rj_abundance[self$traits$index_matrix],
                    self$traits$Ra_abundance[self$traits$index_matrix],
                    self$traits$I2j_abundance[self$traits$index_matrix],
                    self$traits$I2a_abundance[self$traits$index_matrix]),
                  nrow = 8, byrow = TRUE)
    carrying_capacity <- round(self$traits$density_max * self$sim$environment$current$habitat_suitability)[self$traits$index_matrix]
    transformed <- dispersal_fun(
      r = 1,
      tm = 1,
      carrying_capacity = carrying_capacity,
      segment_abundance = segment_abundance
    )
    self$traits$Sj_abundance[self$traits$index_matrix] <- transformed[1, ]
    self$traits$Sa_abundance[self$traits$index_matrix] <- transformed[2, ]
    self$traits$I1j_abundance[self$traits$index_matrix] <- transformed[3, ]
    self$traits$I1a_abundance[self$traits$index_matrix] <- transformed[4, ]
    self$traits$Rj_abundance[self$traits$index_matrix] <- transformed[5, ]
    self$traits$Ra_abundance[self$traits$index_matrix] <- transformed[6, ]
    self$traits$I2j_abundance[self$traits$index_matrix] <- transformed[7, ]
    self$traits$I2a_abundance[self$traits$index_matrix] <- transformed[8, ]
  },
  execution_priority = 3
)

#### Add transition process ####
sim$add_process(
  species = "house_finch",
  process_name = "transition",
  process_fun = function() {
    self$traits$Sa_abundance <- self$traits$Sa_abundance + self$traits$Sj_abundance
    self$traits$Sj_abundance <- matrix(rep(0), nrow = 106, ncol = 161)
    self$traits$I1a_abundance <- self$traits$I1a_abundance + self$traits$I1j_abundance
    self$traits$I1j_abundance <- matrix(rep(0), nrow = 106, ncol = 161)
    self$traits$Ra_abundance <- self$traits$Ra_abundance + self$traits$Rj_abundance
    self$traits$Rj_abundance <- matrix(rep(0), nrow = 106, ncol = 161)
    self$traits$I2a_abundance <- self$traits$I2a_abundance + self$traits$I2j_abundance
    self$traits$I2j_abundance <- matrix(rep(0), nrow = 106, ncol = 161)
  },
  execution_priority = 6
)

#### Add daily disease simulators ####
sim$add_process(
  species = "house_finch",
  process_name = "breeding_season_dynamics",
  process_fun = function() {
    population_new <- daily_siri_summer(
      self$traits$Sj_abundance,
      self$traits$Sa_abundance,
      self$traits$I1j_abundance,
      self$traits$I1a_abundance, 
      self$traits$Rj_abundance, 
      self$traits$Ra_abundance, 
      self$traits$I2j_abundance, 
      self$traits$I2a_abundance,
      self$traits$fecundity,
      self$traits$beta_Sj_winter, 
      self$traits$beta_Sa_winter,
      self$traits$beta_Rj_winter, 
      self$traits$beta_Ra_winter,
      self$traits$recovery_I1j_summer,
      self$traits$recovery_I1a_summer,
      self$traits$recovery_I2j_summer,
      self$traits$recovery_I2a_summer,
      self$traits$mortality_Sj_summer,
      self$traits$mortality_Sa_summer,
      self$traits$mortality_I1j_summer,
      self$traits$mortality_I1a_summer,
      self$traits$mortality_Rj_summer,
      self$traits$mortality_Ra_summer,
      self$traits$mortality_I2j_summer,
      self$traits$mortality_I2a_summer,
      self$sim$environment$current$breeding_season_length,
      self$traits$abundance_threshold,
      self$traits$density_max,
      self$sim$environment$current$habitat_suitability
    )
    if (!exists("population_new") || sum(population_new) == 0) {
      sim$exit()
    } else {
      self$traits$Sj_abundance <- matrix(population_new[1, ], nrow = 106, ncol = 161)
      self$traits$Sa_abundance <- matrix(population_new[2, ], nrow = 106, ncol = 161)
      self$traits$I1j_abundance <- matrix(population_new[3, ], nrow = 106, ncol = 161)
      self$traits$I1a_abundance <- matrix(population_new[4, ], nrow = 106, ncol = 161)
      self$traits$Rj_abundance <- matrix(population_new[5, ], nrow = 106, ncol = 161)
      self$traits$Ra_abundance <- matrix(population_new[6, ], nrow = 106, ncol = 161)
      self$traits$I2j_abundance <- matrix(population_new[7, ], nrow = 106, ncol = 161)
      self$traits$I2a_abundance <- matrix(population_new[8, ], nrow = 106, ncol = 161)
    }
  },
  execution_priority = 7
)

sim$add_process(
  species = "house_finch",
  process_name = "non_breeding_season_dynamics",
  process_fun = function() {
    population_new <- daily_siri_winter(
      self$traits$Sj_abundance,
      self$traits$Sa_abundance,
      self$traits$I1j_abundance,
      self$traits$I1a_abundance, 
      self$traits$Rj_abundance, 
      self$traits$Ra_abundance, 
      self$traits$I2j_abundance, 
      self$traits$I2a_abundance,
      self$traits$beta_Sj_winter, 
      self$traits$beta_Sa_winter,
      self$traits$beta_Rj_winter, 
      self$traits$beta_Ra_winter,
      self$traits$recovery_I1j_winter,
      self$traits$recovery_I1a_winter,
      self$traits$recovery_I2j_winter,
      self$traits$recovery_I2a_winter,
      self$traits$mortality_Sj_winter,
      self$traits$mortality_Sa_winter,
      self$traits$mortality_I1j_winter,
      self$traits$mortality_I1a_winter,
      self$traits$mortality_Rj_winter,
      self$traits$mortality_Ra_winter,
      self$traits$mortality_I2j_winter,
      self$traits$mortality_I2a_winter,
      self$sim$environment$current$breeding_season_length,
      self$traits$abundance_threshold,
      self$traits$density_max,
      self$sim$environment$current$habitat_suitability
    )
    if (sum(population_new) == 0) {
      sim$exit()
    } else {
      self$traits$Sj_abundance <- matrix(population_new[1, ], nrow = 106, ncol = 161)
      self$traits$Sa_abundance <- matrix(population_new[2, ], nrow = 106, ncol = 161)
      self$traits$I1j_abundance <- matrix(population_new[3, ], nrow = 106, ncol = 161)
      self$traits$I1a_abundance <- matrix(population_new[4, ], nrow = 106, ncol = 161)
      self$traits$Rj_abundance <- matrix(population_new[5, ], nrow = 106, ncol = 161)
      self$traits$Ra_abundance <- matrix(population_new[6, ], nrow = 106, ncol = 161)
      self$traits$I2j_abundance <- matrix(population_new[7, ], nrow = 106, ncol = 161)
      self$traits$I2a_abundance <- matrix(population_new[8, ], nrow = 106, ncol = 161)
    }
  },
  execution_priority = 4
)

#### Simulate ####
sim_manager <- metaRangeParallel$new(
  simulation_template = sim,
  generators = list(juvenile_dispersal_gen,
                    adult_dispersal_gen,
                    abundance_gen),
  sample_data = sample_data,
  parallel_threads = 22,
  results_dir = results_dir,
  seed = random_seed,
  species_name = "house_finch"
)

tic()
sim_log <- sim_manager$run()
toc()
