library(raster)
library(terra)
library(tidyverse)
library(poems)
library(epizootic)
library(qs)
library(data.table)
data_dir <- "/Users/caryinstitute/Documents/mgsim/Data/Input"
parallel_cores <- 1
nsims <- 1
burn_in_steps <- 2
timesteps <- 54 + burn_in_steps
random_seed <- 72
results_dir <- "~/Documents/mgsim/Data/Output/epizootic_test"
region <- data_dir %>% file.path("finch_region.qs") %>% qread()
env_corr <- SpatialCorrelation$new(region = region,
                                   amplitude = 0.99,
                                   breadth = 850,
                                   distance_scale = 1000)
env_corr$calculate_compact_decomposition(decimals = 4)
burn_in <- stack(replicate(burn_in_steps, raster(
    file.path(data_dir, "breeding_season_length.tif"), layer = 1
  )))
bsl_raster <- data_dir %>%
  file.path("breeding_season_length.tif") %>%
  stack(bands = 1:54) %>%
  stack(burn_in, .) |>
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
  transmission = rep(0, 4), # because this is pre-outbreak
  transmission_unit = rep(0, 4),
  transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0),
  recovery = rep(0, 4), # pre-outbreak
  recovery_unit = rep(0, 4),
  recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1),
  breeding_season_length = bsl_raster,
  simulation_order = list(c("transition", "season_functions", "results"),
                          c("dispersal", "season_functions", "results")),
  results_selection = c("abundance"),
  results_breakdown = "stages",
  season_functions = list(siri_model_summer, siri_model_winter),
  hs_file = "habitat_suitability_v1",
  mask_file = "native_range_mask",
  verbose = FALSE,
  # I will need to modify this for the next simulation phase
  attribute_aliases = list(
    mortality_Sj_summer = "mortality$summer$a",
    mortality_Sa_summer = "mortality$summer$b",
    mortality_I1j_summer = "mortality$summer$c",
    mortality_I1a_summer = "mortality$summer$d",
    mortality_Rj_summer = "mortality$summer$e",
    mortality_Ra_summer = "mortality$summer$f",
    mortality_I2j_summer = "mortality$summer$g",
    mortality_I2a_summer = "mortality$summer$h",
    mortality_Sj_winter = "mortality$winter$a",
    mortality_Sa_winter = "mortality$winter$b",
    mortality_I1j_winter = "mortality$winter$c",
    mortality_I1a_winter = "mortality$winter$d",
    mortality_Rj_winter = "mortality$winter$e",
    mortality_Ra_winter = "mortality$winter$f",
    mortality_I2j_winter = "mortality$winter$g",
    mortality_I2a_winter = "mortality$winter$h",
    dispersal1 = "dispersal$a",
    dispersal2 = "dispersal$b"
  )
)
capacity_gen <- Generator$new(
  description = "capacity",
  region = region,
  generate_rasters = FALSE,
  # use but don't generate
  burn_in_steps = burn_in_steps,
  time_steps = timesteps,
  generative_requirements = list(
    hs_raster = "file",
    mask_raster = "file",
    initial_abundance = "function",
    carrying_capacity = "function"
  ),
  inputs = c("density_max", "hs_file", "mask_file", "initial_release"),
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
  param = "mask_raster",
  path_template = file.path(data_dir, "%s.tif"),
  path_params = "mask_file",
  file_type = "tif"
)
# Here we subset the hs_matrix to have only the region cells, and we add the burn in
# Also, we tell the generator to generate the carrying_capacity based on "density_max" and "hs_matrix".
capacity_gen$add_function_template(
  param = "carrying_capacity",
  function_def = function(params) {
    hs_matrix <- params$hs_raster %>% as.matrix() %>%
      .[params$region$region_indices, 1:(params$time_steps - params$burn_in_steps)]
    hs_matrix[!is.finite(hs_matrix)] <- 0
    # repeat the first timestep n times as burn in
    hs_matrix <- cbind(replicate(params$burn_in_steps, hs_matrix[, 1]), hs_matrix)
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
    hs_matrix <- params$hs_raster %>%
                 raster::mask(params$mask_raster) %>%
                 as.matrix() %>%
                 .[params$region$region_indices, ] %>%
                 identity()
    hs_matrix[!is.finite(hs_matrix)] <- 0
    hs_vector <- hs_matrix[, 1]
    hs_vector <- round(params$density_max * hs_vector)
    hs_vector[3009] <- params$initial_release
    initial_matrix <- matrix(0, nrow = 8, ncol = 6355)
    initial_matrix[2, ] <- hs_vector
    return(initial_matrix)
  },
  call_params = c("hs_raster", "mask_raster", "density_max", "initial_release", "region")
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
                 dispersal_juv = "dispersal_data")
)
juvenile_dispersal_gen$distance_data <- distance_data
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

# # Transmission parameters
# lhs_generator$set_uniform_parameter("beta_Sa_winter", lower = 0, upper = 0.07588)
# lhs_generator$set_uniform_parameter("beta_Sa_summer", lower = 0, upper = 0.007784)
# lhs_generator$set_triangular_parameter("Sj_multiplier", lower = 0, upper = 8.5,
#                                        mode = 3)
# lhs_generator$set_beta_parameter("beta_I2_modifier", alpha = 1.547023,
#                                  beta = 0.4239236)

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

# # Recovery parameters
# lhs_generator$set_beta_parameter("recovery_I1", alpha = 9.347533,
#                                  beta = 620.1732)
# lhs_generator$set_beta_parameter("recovery_I2", alpha = 1.181112,
#                                  beta = 29.18489)

sample_data <- lhs_generator$generate_samples(number = nsims,
                                              random_seed = random_seed) |>
    mutate(sample = 1:nsims, mortality_Sa_summer = 0,
           mortality_I2j_summer = mortality_I2_modifier * mortality_I1j_summer,
           mortality_I2j_winter = mortality_I2_modifier * mortality_I1j_winter,
           mortality_I2a_winter = mortality_I2_modifier * mortality_I1a_winter,
           mortality_I2a_summer = mortality_I2_modifier * mortality_I1a_summer,
           mortality_Rj_summer = mortality_Sj_summer,
           mortality_Ra_summer = mortality_Sa_summer,
           mortality_Rj_winter = mortality_Sj_winter,
           mortality_Ra_winter = mortality_Sa_winter)
handler <- SimulationHandler$new(
  sample_data = sample_data,
  model_template = model_template,
  generators = list(juvenile_dispersal_gen,
                    adult_dispersal_gen,
                    capacity_gen),
  parallel_cores = parallel_cores,
  results_dir = results_dir
)
self <- handler
j <- 1
# Check for error messages
if (!is.null(self$error_messages)) {
  error_messages <- self$error_messages
  self$error_messages <- NULL
  stop(error_messages, call. = FALSE)
}

# Check that model and sample data is present
if (is.null(self$model_template) | length(self$sample_data) == 0) {
  cli_abort(c("Simulations cannot run unless there is a `model_template`
                    or a `sample_data` with at least one row."))
}

# Check that model simulator/function is present
if (is.null(self$model_simulator)) {
  stop("The model simulator has not been set", call. = FALSE)
} else if (is.null(self$model_simulator$simulation_function)) {
  stop("The model simulator function has not been set", call. = FALSE)
}

# Check that the results directory is present and exists
if (is.null(self$results_dir)) {
  stop("No output directory set for results", call. = FALSE)
}
if (!dir.exists(self$results_dir)) {
  stop(paste("Could not find results directory", self$results_dir), call. = FALSE)
}
if (is.null(self$results_ext)) {
  self$results_ext <- ".qs" # reinstate default
}

# Create a nested simulation (or descendant) model for cloning
self$nested_model <- self$model_template$new_clone(template = self$model_template)

# Allow extra attachments to be passed
if ("nested_model" %in% names(self$attached)) {
  self$nested_model$attached <- self$attached$nested_model
}

model_sample_columns <- which(names(self$sample_data) %in% self$nested_model$get_attribute_aliases())
if (length(model_sample_columns) > 0) {
  self$nested_model$attached$sample_model_names <- names(self$sample_data)[model_sample_columns]
  self$nested_model$sample_attributes <- self$nested_model$attached$sample_model_names
}
if (!is.null(self$generators)) {
  self$nested_model$attached$sample_generative_names <- list()
  dispersal_count <- 0

  for (i in 1:length(self$generators)) {
    generator <- self$generators[[i]]

    if ("DispersalGenerator" %in% class(generator)) {
      dispersal_count <- dispersal_count + 1
      self$nested_model$attached$sample_generative_names[[i]] <- paste0("dispersal", dispersal_count)
      self$nested_model$sample_attributes <- unique(c(self$nested_model$sample_attributes, paste0("dispersal", dispersal_count)))
    } else {
      self$nested_model$attached$sample_generative_names[[i]] <- generator$outputs
      self$nested_model$sample_attributes <- unique(c(self$nested_model$sample_attributes, generator$outputs))
    }
  }
}

# Check the completeness/consistency of the first sample only
model <- self$nested_model$clone()
self$set_model_sample(model, 1)
inputs <- model
# Check that all inputs are valid
inputs <- check_simulator_inputs(inputs)
list2env(inputs, envir = environment())

population_abundance <- colSums(initial_abundance)

# Simulator reference object for dynamic attachments and results
# (accessed via user-defined functions)
simulator <- SimulatorReference$new()

# Generate simulation functions
transition_function <- disease_transitions(stages, compartments)

translocation_function <- population_transformation(
  replicates = replicates,
  time_steps = time_steps,
  populations = populations,
  demographic_stochasticity = demographic_stochasticity,
  density_stages = density_stages,
  transformation = inputs[["translocation"]],
  simulator = simulator,
  name = "translocation"
)
harvest_function <- population_transformation(
  replicates = replicates,
  time_steps = time_steps,
  populations = populations,
  demographic_stochasticity = demographic_stochasticity,
  density_stages = density_stages,
  transformation = inputs[["harvest"]],
  simulator = simulator,
  name = "harvest"
)

if ("mortality" %in% simulation_order) {
  mortality_function <- population_transformation(
    replicates = replicates,
    time_steps = time_steps,
    populations = populations,
    demographic_stochasticity = demographic_stochasticity,
    density_stages = density_stages,
    transformation = inputs[["mortality_function"]],
    simulator = simulator,
    name = "mortality"
  )
}

if (length(inputs[["dispersal"]]) == 1) {
  args <- list(
    "replicates" = replicates,
    "time_steps" = time_steps,
    "years_per_step" = 1,
    "populations" = populations,
    "demographic_stochasticity" = demographic_stochasticity,
    "density_stages" = density_stages,
    "dispersal" = inputs[["dispersal"]],
    "dispersal_stages" = dispersal_stages,
    "simulator" = simulator
  )

  # Check if each argument exists and add it to the args list if it does
  if (exists("dispersal_source_n_k")) args$dispersal_source_n_k <- dispersal_source_n_k
  if (exists("dispersal_target_k")) args$dispersal_target_k <- dispersal_target_k
  if (exists("dispersal_target_n")) args$dispersal_target_n <- dispersal_target_n
  if (exists("dispersal_target_n_k")) args$dispersal_target_n_k <- dispersal_target_n_k

  # Call the function with the args list
  do.call(population_dispersal, args)
} else if (length(inputs[["dispersal"]]) > 1) {
  dispersal_functions = list()
  for (d in 1:length(inputs[["dispersal"]])) {
    args <- list(
      "replicates" = replicates,
      "time_steps" = time_steps,
      "years_per_step" = 1,
      "populations" = populations,
      "demographic_stochasticity" = demographic_stochasticity,
      "density_stages" = density_stages,
      "dispersal" = inputs[["dispersal"]][[d]],
      "dispersal_stages" = dispersal_stages,
      "simulator" = simulator
    )

    # Check if each argument exists and add it to the args list if it does
    if (exists("dispersal_source_n_k")) args$dispersal_source_n_k <- dispersal_source_n_k
    if (exists("dispersal_target_k")) args$dispersal_target_k <- dispersal_target_k
    if (exists("dispersal_target_n")) args$dispersal_target_n <- dispersal_target_n
    if (exists("dispersal_target_n_k")) args$dispersal_target_n_k <- dispersal_target_n_k

    # Call the function with the args list
    dispersal_functions[[d]] <- do.call(population_dispersal, args)
  }
}

if (exists("season_functions")) {
  season_function_list <- list()
  for (i in 1:length(season_functions)) {
    season_function_list[[i]] <- disease_transformation(
      list(
        "replicates" = replicates,
        "time_steps" = time_steps,
        "seasons" = seasons,
        "compartments" = compartments,
        "populations" = populations,
        "demographic_stochasticity" = demographic_stochasticity,
        "stages" = stages,
        "abundance_threshold" = abundance_threshold,
        "mortality" = mortality[[i]],
        "mortality_unit" = mortality_unit[[i]],
        "fecundity" = fecundity[[i]],
        "fecundity_unit" = fecundity_unit[[i]],
        "fecundity_mask" = fecundity_mask[[i]],
        "transmission" = transmission[[i]],
        "transmission_unit" = transmission_unit[[i]],
        "transmission_mask" = transmission_mask[[i]],
        "recovery" = recovery[[i]],
        "recovery_unit" = recovery_unit[[i]],
        "recovery_mask" = recovery_mask[[i]],
        "transformation" = season_functions[[i]],
        "simulator" = simulator,
        "name" = paste0("season", i, collapse = "_")
      )
    )
  }
}

result_functions <- disease_results(
  replicates = replicates,
  time_steps = time_steps,
  seasons = seasons,
  stages = stages,
  compartments = compartments,
  coordinates = inputs[["coordinates"]],
  initial_abundance = initial_abundance,
  results_selection = results_selection,
  results_breakdown = results_breakdown
)
results_list <- result_functions$initialize_attributes()

r <- 1
tm <- 1
season <- 2

segment_abundance <- initial_abundance
population_abundance <- .colSums(initial_abundance,
                                 m = segments, n = populations)
occupied_indices <- which(as.logical(population_abundance))
occupied_populations <- length(occupied_indices)

indices <- map(1:stages, \(s) seq(s, segments, stages))
stage_abundance <- segment_abundance[indices[[s]],]

list2env(args, environment())

if (is.list(dispersal) && is.data.frame(dispersal[[1]]) && nrow(dispersal[[1]])) { # compact matrix data in list (as per DispersalModel class)

  # Unpack dispersal data and determine compact matrix dimensions
  dispersal_data <- dispersal[[1]]
  dispersal_compact_rows <- max(dispersal_data[, c("emigrant_row", "immigrant_row")])

  # Are dispersals changing over time?
  dispersals_change_over_time <- (length(dispersal) > 1)
  if (dispersals_change_over_time) {
    dispersal_data_changes <- dispersal
    dispersal_data_changes[[1]] <- dispersal_data_changes[[1]][NULL,]
  }

}

# Release dispersal from memory
dispersal <- NULL

# Create a compact matrix of dispersal rates
dispersal_compact_matrix <- array(0, c(dispersal_compact_rows, populations))
dispersal_compact_matrix[as.matrix(dispersal_data[, c("emigrant_row", "source_pop")])] <- dispersal_data$dispersal_rate

# Does dispersal depend on source population abundance N divided by carrying capacity K?
dispersal_depends_on_source_pop_n_k <- (is.list(dispersal_source_n_k) && (is.numeric(dispersal_source_n_k$cutoff) ||
                                                                            is.numeric(dispersal_source_n_k$threshold)))
dispersal_depends_on_target_pop_n_k <- (is.list(dispersal_target_n_k) && (is.numeric(dispersal_target_n_k$threshold) ||
                                                                            is.numeric(dispersal_target_n_k$cutoff)))
# Setup density dependence dispersal parameters
if (dispersal_depends_on_source_pop_n_k) {

  # Convert NULL to zero in source N/K cutoff or one in threshold
  if (dispersal_depends_on_source_pop_n_k) {
    if (is.null(dispersal_source_n_k$cutoff)) dispersal_source_n_k$cutoff <- 0
    if (is.null(dispersal_source_n_k$threshold)) dispersal_source_n_k$threshold <- 1
  }

  # Check threshold > cutoff
  if (dispersal_source_n_k$threshold <= dispersal_source_n_k$cutoff) {
    dispersal_depends_on_source_pop_n_k <- FALSE
    warning("Dispersal density dependence for source N/K threshold must be greater than cutoff => not used", call. = FALSE)
  }
}
  if (dispersal_depends_on_target_pop_n_k) {

    # Convert NULL to zero in target N/K threshold or cutoff
    if (is.null(dispersal_target_n_k$threshold)) dispersal_target_n_k$threshold <- 0
    if (is.null(dispersal_target_n_k$cutoff)) dispersal_target_n_k$cutoff <- 0
  }

  # Create a map of compact array indices for mapping dispersers (emigrants) to target populations
  dispersal_target_pop_map <- array(0, c(dispersal_compact_rows, populations))
  dispersal_target_pop_map[as.matrix(dispersal_data[, c("emigrant_row", "source_pop")])] <- dispersal_data$target_pop


dispersal_compact_indices <- array(1:(dispersal_compact_rows*populations), c(dispersal_compact_rows, populations))
dispersal_immigrant_map <- array(0, c(dispersal_compact_rows, populations))
dispersal_immigrant_map[as.matrix(dispersal_data[, c("emigrant_row", "source_pop")])] <- dispersal_compact_indices[as.matrix(dispersal_data[, c("immigrant_row", "target_pop")])]

# Release variables from memory
dispersal_data <- NULL; dispersal_compact_indices <- NULL

# Calculate occupied population number
occupied_populations <- length(occupied_indices)

# Apply any spatio-temporal dispersal changes
dispersal_compact_matrix_tm <- simulator$attached$dispersal_compact_matrix_tm
if (tm == 1 || !dispersals_change_over_time) {
  dispersal_compact_matrix_tm <- dispersal_compact_matrix
} else if (dispersals_change_over_time && nrow(dispersal_data_changes[[tm]])) { # and tm > 1
  dispersal_compact_matrix_tm[as.matrix(dispersal_data_changes[[tm]][, c("emigrant_row","source_pop")])] <- dispersal_data_changes[[tm]]$dispersal_rate
}
simulator$attached$dispersal_compact_matrix_tm <- dispersal_compact_matrix_tm

# Select dispersals for occupied populations
occupied_dispersals <- dispersal_compact_matrix_tm[, occupied_indices]

# Calculate density abundance
if (dispersal_depends_on_source_pop_n_k || dispersal_depends_on_target_pop_n_k) {
  density_abundance <- .colSums(stage_abundance*as.numeric(density_stages), m = length(density_stages), n = populations)
}

# Modify dispersal rates when dispersal depends on source population N/K
if (dispersal_depends_on_source_pop_n_k) {

  # Density dependent multipliers
  dd_multipliers <- array(1, populations)

  # Calculate the source N/K multipliers
  abundance_on_capacity <- density_abundance/carrying_capacity
  dd_multipliers[which(abundance_on_capacity <= dispersal_source_n_k$cutoff)] <- 0
  modify_pop_indices <- which(carrying_capacity > 0 & dd_multipliers > 0 &
                                abundance_on_capacity < dispersal_source_n_k$threshold)
  dd_multipliers[modify_pop_indices] <- ((abundance_on_capacity[modify_pop_indices] -
                                            array(dispersal_source_n_k$cutoff, populations)[modify_pop_indices])/
                                           array(dispersal_source_n_k$threshold - dispersal_source_n_k$cutoff,
                                                 populations)[modify_pop_indices]*
                                           dd_multipliers[modify_pop_indices])

  # Apply modifying multipliers to dispersals
  occupied_dispersals <- (occupied_dispersals*matrix(dd_multipliers[occupied_indices],
                                                     nrow = dispersal_compact_rows,
                                                     ncol = occupied_populations, byrow = TRUE))

} # dispersal depends on source pop N/K?

# Select occupied dispersal non-zero indices
occupied_dispersal_indices <- which(as.logical(occupied_dispersals)) # > 0

# Modify dispersal rates when dispersal depends on target population K, N, or N/K
if (dispersal_depends_on_target_pop_n_k) {

  # Density dependent multipliers
  dd_multipliers <- array(1, populations)

  # Calculate the target N/K multipliers
  if (dispersal_depends_on_target_pop_n_k) {
    dd_multipliers[which(carrying_capacity <= 0)] <- 0
    abundance_on_capacity <- density_abundance/carrying_capacity
    if (all(dispersal_target_n_k$threshold < dispersal_target_n_k$cutoff)) { # overcrowded cell avoidance \
      dd_multipliers[which(abundance_on_capacity >= dispersal_target_n_k$cutoff)] <- 0
      modify_pop_indices <- which(abundance_on_capacity > dispersal_target_n_k$threshold & dd_multipliers > 0)
      dd_multipliers[modify_pop_indices] <- ((array(dispersal_target_n_k$cutoff, populations)[modify_pop_indices] -
                                                abundance_on_capacity[modify_pop_indices])/
                                               array(dispersal_target_n_k$cutoff - dispersal_target_n_k$threshold,
                                                     populations)[modify_pop_indices]*
                                               dd_multipliers[modify_pop_indices])
    } else if (all(dispersal_target_n_k$threshold > dispersal_target_n_k$cutoff)) { # seek company /
      dd_multipliers[which(abundance_on_capacity <= dispersal_target_n_k$cutoff)] <- 0
      modify_pop_indices <- which(abundance_on_capacity < dispersal_target_n_k$threshold & dd_multipliers > 0)
      dd_multipliers[modify_pop_indices] <- ((abundance_on_capacity[modify_pop_indices] -
                                                array(dispersal_target_n_k$cutoff, populations)[modify_pop_indices])/
                                               array(dispersal_target_n_k$threshold - dispersal_target_n_k$cutoff,
                                                     populations)[modify_pop_indices]*
                                               dd_multipliers[modify_pop_indices])
    }
  }

  # Select multipliers via target populations for non-zero occupied dispersals
  selected_dd_multipliers <- dd_multipliers[dispersal_target_pop_map[, occupied_indices][occupied_dispersal_indices]]

  # Apply modifying multipliers to dispersals
  modify_indices <- which(selected_dd_multipliers < 1)
  if (length(modify_indices)) {
    modify_dipersal_indices <- occupied_dispersal_indices[modify_indices]
    occupied_dispersals[modify_dipersal_indices] <- occupied_dispersals[modify_dipersal_indices]*selected_dd_multipliers[modify_indices]
    occupied_dispersal_indices <- which(as.logical(occupied_dispersals)) # > 0
  }

} # dispersal depends on target pop N, K or N/K?

# Perform dispersal for each participating stage
for (stage in which(dispersal_stages > 0)) {

  # Disperser generation via abundance and corresponding dispersal rates
  occupied_abundance <- stage_abundance[stage, occupied_indices]
  occupied_abundance_rep <- stage_abundance[rep(stage, dispersal_compact_rows), occupied_indices]
  dispersers <- array(0, c(dispersal_compact_rows, occupied_populations))

  # Generate dispersers
  if (demographic_stochasticity) { # via binomial distribution
    dispersers[occupied_dispersal_indices] <- stats::rbinom(length(occupied_dispersal_indices), occupied_abundance_rep[occupied_dispersal_indices],
                                                            occupied_dispersals[occupied_dispersal_indices]*dispersal_stages[stage])
  } else { # deterministic
    dispersers[occupied_dispersal_indices] <- round(occupied_abundance_rep[occupied_dispersal_indices]*
                                                      occupied_dispersals[occupied_dispersal_indices]*dispersal_stages[stage])
  }

  # Calculate emigrants
  emigrants <- array(.colSums(dispersers, m = dispersal_compact_rows, n = occupied_populations))

  # Check consistency of emigrants (not to exceed abundances)
  excessive_indices <- which(emigrants > occupied_abundance)
  if (length(excessive_indices) > 0) { # reduce emigrants to equal abundance via random sampling
    for (excessive_index in excessive_indices) {
      excessive_rows <- which(as.logical(dispersers[, excessive_index])) # > 0
      excessive_dispersers <- dispersers[excessive_rows, excessive_index]
      disperser_reduction <- emigrants[excessive_index] - occupied_abundance[excessive_index]
      for (remove_row_index in rep(excessive_rows,
                                   times = excessive_dispersers)[sample(sum(excessive_dispersers),
                                                                        size = disperser_reduction)]) {
        dispersers[remove_row_index, excessive_index] <- dispersers[remove_row_index, excessive_index] - 1
      }
    }
    emigrants[excessive_indices] <- occupied_abundance[excessive_indices]
  }

  # Update occupied stage abundance
  stage_abundance[stage, occupied_indices] <- stage_abundance[stage, occupied_indices] - emigrants

  # Calculate immigrants via dispersal immigrant map
  disperser_indices <- which(as.logical(dispersers)) # > 0
  immigrant_array <- array(0, c(dispersal_compact_rows, populations))
  immigrant_array[dispersal_immigrant_map[, occupied_indices][disperser_indices]] <- dispersers[disperser_indices]
  immigrants <- .colSums(immigrant_array, m = dispersal_compact_rows, n = populations)

  # Update population abundances
  stage_abundance[stage,] <- stage_abundance[stage,] + immigrants

}

# Perform additional dispersal for overcrowded cells (only to cells with room)
if ((dispersal_depends_on_target_pop_n_k && all(dispersal_target_n_k$threshold < dispersal_target_n_k$cutoff))) {

  depends_on_target_pop_n_k <- (dispersal_depends_on_target_pop_n_k && all(dispersal_target_n_k$threshold < dispersal_target_n_k$cutoff))

  # Get all updated dispersal rates
  dispersals <- dispersal_compact_matrix_tm
  dispersals[, occupied_indices] <- occupied_dispersals

  # Identify overcrowded cells based on stages affected by density
  density_abundance <- .colSums(stage_abundance*as.numeric(density_stages), m = length(density_stages), n = populations)
  stage_indices <- which(density_stages > 0 & dispersal_stages > 0)
  excessive_indices <- c()

  if (depends_on_target_pop_n_k) {
    excessive_indices <- unique(c(excessive_indices,
                                  which(density_abundance/carrying_capacity > dispersal_target_n_k$cutoff)))
  }

  # Disperse excess from each overcrowded cell (in random order)
  for (excessive_index in excessive_indices[sample(length(excessive_indices))]) {

    # Determine dispersal targets and rates with enough room for excess from overcrowded cell
    dispersal_indices <- which(dispersals[, excessive_index] > 0)
    target_indices <- dispersal_target_pop_map[, excessive_index][dispersal_indices]
    if (depends_on_target_pop_n_k) {
      indices_with_room <- which(((density_abundance + 1)/carrying_capacity <= dispersal_target_n_k$cutoff)[target_indices])
    }
    dispersal_indices <- dispersal_indices[indices_with_room]
    target_indices <- target_indices[indices_with_room]

    # Disperse excess one at a time sampled via the cell stage abundance distribution
    rep_stage_indices <- rep(stage_indices, times = stage_abundance[stage_indices, excessive_index])
    abundance_excess <- 0
    if (depends_on_target_pop_n_k) {
      abundance_excess <- max(abundance_excess, density_abundance[excessive_index] - floor(dispersal_target_n_k$cutoff*carrying_capacity[excessive_index]))
    }
    for (stage_i in rep_stage_indices[sample(1:length(rep_stage_indices), size = abundance_excess)]) {
      if (length(target_indices)) {

        # Sample target cell
        target_i <- target_indices[sample(length(target_indices), size = 1,
                                          prob = dispersals[dispersal_indices, excessive_index])]

        # Perform dispersal
        stage_abundance[stage_i, excessive_index] <- stage_abundance[stage_i, excessive_index] - 1 # emigrant
        stage_abundance[stage_i, target_i] <- stage_abundance[stage_i, target_i] + 1 # immigrant

        # Update target density abundance and potential targets if it becomes full
        density_abundance[target_i] <- density_abundance[target_i] + 1
        if ((depends_on_target_pop_n_k && density_abundance[target_i]/carrying_capacity[target_i] >= dispersal_target_n_k$cutoff)) { # remove from potential targets
          full_index <- which(target_indices == target_i)
          target_indices <- target_indices[-full_index]
          dispersal_indices <- dispersal_indices[-full_index]
        }
      }
    }
  }
}
