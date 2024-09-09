self <- handler
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
if (length(model$error_messages)) {
  stop(c("Error(s) setting model sample attributes: ", model$error_messages), call. = FALSE)
}
if (!model$is_complete()) {
  incomplete_message <- "Model attributes are incomplete"
  if (!model$is_consistent()) {
    incomplete_message <- paste(incomplete_message, "/inconsistent", sep = "")
  }
  incomplete_message <- paste0(incomplete_message, ": ", paste(model$incomplete_attributes(), collapse = ", "))
  stop(incomplete_message, call. = FALSE)
}
model <- NULL # release from memory
# Clone the model
model <- self$nested_model$clone()

# Set the model sample attributes
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

if ("mortality" %in% flatten(simulation_order)) {
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

if ("dispersal" %in% flatten(simulation_order)) {
  dispersal_function <- disease_dispersal(
    replicates,
    time_steps,
    populations,
    demographic_stochasticity,
    dispersal,
    dispersal_type,
    dispersal_source_n_k = if (exists("dispersal_source_n_k")) dispersal_source_n_k else NULL,
    dispersal_target_k = if (exists("dispersal_target_k")) dispersal_target_k else NULL,
    dispersal_target_n = if (exists("dispersal_target_n")) dispersal_target_n else NULL,
    dispersal_target_n_k = if (exists("dispersal_target_n_k")) dispersal_target_n_k else NULL,
    stages = stages,
    compartments = compartments,
    simulator
  )
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
