library(checkmate)
library(furrr)

#' Interpolation of Missing Timesteps in a Raster Stack
#'
#' Interpolates missing time points in a raster stack with layers for different time points. Intended
#' for `terra::SpatRaster` use only.
#'
#' @importFrom terra rast approximate xmin xmax ymin ymax crs res values<-
#' @import purrr
#' @param raster_stack A `terra::SpatRaster` with at least two layers.
#' @param source_time A numeric vector indicating the timesteps of the `raster_stack`, whatever those
#' may be, e.g., `c(1940, 1945)`.
#' @param target_time A numeric vector indicating the series of timesteps desired in the output, e.g.,
#' `c(1940, 1941, 1942. 1943, 1944, 1945)`.
#' @param time_label `terra::SpatRaster`s do not allow raw numbers as names of raster layers. Therefore
#' a string is needed to place before the timestep number, e.g., "BP", "BCE".
#' @param ... Does nothing. A placeholder for future code improvements.
#' @param method Interpolation method passed on to `terra::approximate`. Default is "linear";
#' alternative is "constant", which implements a step function.
#' @return A `terra::SpatRaster` with as many layers as the length of `target_time`.
#' @export
interpolate_raster <- function(raster_stack,
                               source_time,
                               target_time,
                               time_label,
                               ...,
                               method = "linear") {
  if (!inherits(raster_stack, "SpatRaster")) {
    stop("raster_stack must be a terra::SpatRaster")
  }
  if (!all(is.numeric(source_time), is.numeric(target_time))) {
    stop("source time and target time must be numeric vectors")
  }
  if (length(target_time) < length(source_time)) {
    stop(
      "Target time vector must be longer than source time vector
         \n (otherwise there's no interpolation)"
    )
  }
  if (!is.character(time_label)) {
    stop("time_label must be a string")
  }
  template <- raster_stack[[1]]
  # Create output raster stack
  outputStack <- rast(
    nlyrs = length(target_time),
    nrows = nrow(template),
    ncols = ncol(template),
    xmin = xmin(template),
    xmax = xmax(template),
    ymin = ymin(template),
    ymax = ymax(template),
    crs = crs(template),
    resolution = res(template),
    vals = NA,
    names = target_time %>% map_chr(~ paste0(time_label, .))
  )
  for (i in seq_along(source_time)) {
    j <- which(target_time == source_time[i])
    time <- target_time[j]
    inSeq <- which(source_time == time)
    if (length(inSeq) != 0) {
      r <- raster_stack[[inSeq]]
      values(outputStack[[j]]) <- r[]
    }
  }
  # urbanStack <- setZ(urbanStack, z = ts, name = "Years BP")
  interpolated_urban <- approximate(outputStack, method = method)
  return(interpolated_urban)
}

ZeroFillPFW <- function(InputData, SpeciesCodes, rollup = TRUE) {
  #The first three steps collect information that
  # will be necessary for taxonomic roll-up and zero-filling. First, we create
  # a table that contains only the information describing the observation
  # periods (and not anything specific to any single birds species). Second,
  # we need to count the number of bird species for which zero-filling has been
  # requested, in order to know how many times we will need to loop through the
  # zero-filling process. Third, we need to check the species codes for which
  # zero-filling was requested by comparing these to a list of all species codes
  # found in InputData.

  #Create the table of sampling event information from all observation periods
  PFW_SED <- InputData %>%
    select(-OBS_ID,
           -SPECIES_CODE,
           -HOW_MANY,
           -PLUS_CODE,
           -VALID,
           -REVIEWED) %>%
    distinct()

  #Count the number of species codes
  NSpecies <- length(SpeciesCodes)

  #Identify any species codes that are not found in the data
  InvalidSppCodes <- InputData %>%
    select(SPECIES_CODE) %>%
    distinct() %>%
    anti_join(y = .,
              x = tibble(SPECIES_CODE = SpeciesCodes),
              by = "SPECIES_CODE")
  #For each non-existing species code, print an error message
  if (nrow(InvalidSppCodes) >= 1) {
    for (i in 1:nrow(InvalidSppCodes)) {
      cat(
        paste(
          "THE SPECIES CODE \'",
          InvalidSppCodes$SPECIES_CODE[i],
          "\' DOES NOT EXIST\n",
          sep = ""
        )
      )
      cat(" IN THE DATA TABLE, AND ZERO-FILLING IS NOT POSSIBLE.\n\n")
    }
    #break out of the function by returning a null result
    cat("THIS FUNCTION IS STOPPING BECAUSE ZERO-FILLING OF ALL\n")
    cat(" REQUESTED SPECIES CODES IS NOT POSSIBLE.\n\n")
    return(NULL)
  }

  #Loop through the species codes in turn, possibly doing a taxonomic roll-up,
  # and always zero-filling each species' data.
  for (i in 1:NSpecies) {
    #
    #Check if this species code is "nobird" and if so skip over this
    # species code.
    if (SpeciesCodes[i] == "nobird") {
      cat("THE SPECIES CODE \'nobird\' INDICATES AN ABSENCE OF BIRDS,\n")
      cat(" AND IS NOT AN ACTUAL BIRD SPECIES. NO ZERO-FILLING WILL\n")
      cat("BE DONE FOR THE SPECIES_CODE \'nobird\'.\n\n")
      next
    }
    #Place the SPECIES_CODE for the ith species into Spp_to_Zero_Fill.
    Spp_to_Zero_Fill <- SpeciesCodes[i]
    #
    #If a taxonomic roll-up is to be done, then do so for the focal species.
    if (rollup) {
      InputData <- InputData %>%
        mutate(SPECIES_CODE = case_when(
          (is.na(alt_full_spp_code)
           ~ SPECIES_CODE),
          ((
            !is.na(alt_full_spp_code) &
              alt_full_spp_code == Spp_to_Zero_Fill
          )
          ~ alt_full_spp_code
          ),
          ((
            !is.na(alt_full_spp_code) &
              alt_full_spp_code != Spp_to_Zero_Fill
          )
          ~ SPECIES_CODE
          )
        ))
    }

    #Create template for zero-count records for the focal species.
    FakeSppSpecificFields <- as_tibble_row(
      list(
        OBS_ID = "OBSnull",
        SPECIES_CODE = Spp_to_Zero_Fill,
        HOW_MANY = 0,
        PLUS_CODE = NA,
        VALID = 1,
        REVIEWED = 0
      )
    )

    #Zero-fill the data for the focal species.
    ZeroFilledFocalSpecies <- InputData %>%
      filter(SPECIES_CODE == Spp_to_Zero_Fill) %>%
      bind_rows(., bind_cols(PFW_SED, FakeSppSpecificFields)) %>%
      group_by(across(c(names(PFW_SED), "SPECIES_CODE"))) %>%
      summarize(
        OBS_ID = min(OBS_ID),
        HOW_MANY = sum(HOW_MANY),
        PLUS_CODE = max(PLUS_CODE),
        VALID = min(VALID),
        REVIEWED = max(REVIEWED)
      ) %>%
      ungroup()
    ZeroFilledFocalSpecies$PLUS_CODE <- as.logical(ZeroFilledFocalSpecies$PLUS_CODE)

    #Put the zero-filled data into a table that accumulates data from all of
    # the species named in the vector SpeciesCodes. Either create the tibble
    # to be returned if it does not already exist, or append rows to an
    # existing object.
    if ((exists("OutputData")))
      OutputData <- bind_rows(OutputData, ZeroFilledFocalSpecies)
    else
      OutputData <- ZeroFilledFocalSpecies
  }
  #After all species' data are zero-filled, return the zero-filled data to the
  # user. Unless the user of the function directs the output into some object,
  # the output will be printed to the console and not be saved for subsequent use.
  #Also, if a taxonomic roll-up has been done, print a message to remind users
  # this was done.
  if (rollup) {
    cat("\nNOTE: ANY RECORDS OF BIRDS REPORTED TO THE LEVEL OF A SUBSPECIES OR\n")
    cat(" OTHER RECOGNIZABLE FORM HAVE BEEN TREATED AS MEMBERS OF THE FULL\n")
    cat(" SPECIES FOR THE PURPOSE OF ZERO-FILLING, AND 'species_code' VALUES\n")
    cat(" FOR THESE RECOGNIZABLE FORMS HAVE BEEN REPLACE BY THE CODES FOR THE\n")
    cat(" FULL SPECIES IN THE OUTPUT.\n\n")
  }
  return(OutputData)
}


#' Unpack the contents of a closure into the global environment
#'
#' Sometimes for debugging purposes you need to unpack a closure like a
#' suitcase.
#'
#' @param closure The closure to be unpacked.
#' @return Invisibly copies objects into the global environment.
unpack_closure <- function(closure) {
  closure_env <- environment(closure)
  vars <- ls(envir = closure_env)

  for (var in vars) {
    assign(var, get(var, envir = closure_env), envir = .GlobalEnv)
  }
}

## Write simulation numbers to a batch transfer file
write_sftp_transfer <- function(source_dir,
                                dest_dir,
                                sim_numbers,
                                filename) {
  # Assertions
  assert_character(source_dir, min.chars = 1)
  assert_character(dest_dir, min.chars = 1)
  assert_integerish(
    sim_numbers,
    lower = 1,
    any.missing = FALSE,
    min.len = 1
  )
  assert_path_for_output(filename, overwrite = T)

  # Check if destination directories already exist
  destinations <- map_chr(sim_numbers, function(n) {
    file.path(dest_dir, paste0("simulation", n))
  })
  dir_exists <- map_lgl(destinations, dir.exists)

  instructions <- map_chr(sim_numbers[!dir_exists], function(n) {
    paste0("get -r ", source_dir, "/simulation", n, " ", dest_dir)
  }) |>
    paste(collapse = "\n")
  write(instructions, filename)
}

## Create a weighted ensemble mean of selected samples for a given stage +
## compartment
# Output: an array of 106 rows, 161 columns, and 154 time points (1940-2016,
# winter/summer)
ensemble_mean <- function(samples, weights, compartment, data_dir) {
  # Assertions
  assert_integerish(
    samples,
    lower = 1,
    upper = 10000,
    any.missing = FALSE,
    min.len = 1
  )
  assert_numeric(
    weights,
    lower = 0,
    upper = 1,
    any.missing = F,
    len = length(samples)
  )
  assert_string(compartment, min.chars = 2)
  assert_choice(compartment,
                c("Sj", "Sa", "I1j", "I1a", "Rj", "Ra", "I2j", "I2a"))
  assert_directory_exists(data_dir, access = "r")

  # Read in data
  dir_vec <- map_chr(samples, function(n) {
    file.path(data_dir, paste0("simulation", n))
  })
  num_slashes <- str_count(dir_vec[1], "/")
  year_lookup <- data.frame(index = 1:77, Year = 1940:2016)
  file_info <- dir_vec |>
    lapply(list.files, full.names = TRUE) |>
    lapply(
      function(p)
        data.frame(
          path = p,
          index = as.numeric(str_split_i(p, "/", num_slashes + 2) |>
                               str_split_i("_", 1)),
          season = str_extract(p, "summer|winter")
        ) |>
        mutate(season = factor(season, levels = c(
          "winter", "summer"
        ))) |>
        arrange(index, season) |>
        left_join(year_lookup, by = join_by(index)) |>
        filter(str_detect(path, compartment))
    ) |>
    map(fill_missing_years,
        start_year = 1940,
        end_year = 2016)

  plan(multisession)
  # Read in .qs files and create a weighted ensemble mean
  ensemble <- future_map2(file_info, weights, function(info, weight) {
    # Preallocate a 3D array for sim_array results
    sim_array_results <- array(0, dim = c(106, 161, nrow(info)))
    for (i in seq_len(nrow(info))) {
      # Read the .qs file
      sim_data <- read_or_zero(info$path[i])
      # Fill the 3D array with the data
      sim_array_results[, , i] <- sim_data
    }
    return(sim_array_results * weight)
  })

  weighted_ensemble <- Reduce(`+`, ensemble) / sum(weights)
  return(weighted_ensemble)
}

## Plot priors and posteriors from Approximate Bayesian Computation
plot_abc_posteriors <- function(abc, priors) {
  # Assertions
  assert_class(abc, "abc")
  assert_choice(abc$method, c("rejection", "loclinear", "neuralnet"))
  assert_data_frame(priors, min.rows = 1, min.cols = 1, col.names = "unique")
  assert_true(all.equal(names(priors), colnames(abc$adj.values)))

  prior_samples <- priors |>
    mutate(Distribution = "Prior")
  param_names <- colnames(abc$adj.values)
  param_labels <- c(
    "Proportion of Juveniles Dispersing",
    "Proportion of Adults Dispersing",
    "Maximum Dispersal Distance (Juv.)",
    "Maximum Dispersal Distance (Ad.)",
    "Density over this value\n reduces outward dispersal",
    "Density over this value\n ends outward dispersal",
    "Density over this value\n reduces inward dispersal",
    "Density over this value\n ends inward dispersal",
    "Quasi-extinction threshold",
    "# Finches released in NY",
    "Maximum Population Density (2151 sq. km.)",
    "Fecundity (total breeding season)",
    "# Finches initially infected",
    "Beta, Susceptible Adults, Winter",
    "Beta, Susceptible Adults, Summer",
    "How much is beta reduced\n in subsequent infections",
    "Mortality, Susceptible Juveniles, Winter",
    "Mortality, Susceptible Adults, Winter",
    "Mortality, Susceptible Juveniles, Summer",
    "Mortality, 1st infection Juveniles, Summer",
    "Mortality, 1st infection Juveniles, Winter",
    "Mortality, 1st infection Adults, Summer",
    "Mortality, 1st infection Adults, Winter",
    "How much is mortality reduced\n after 1st infection",
    "Mortality, Susceptible Adults, Summer",
    "Mortality, 2nd+ infection Juveniles, Summer",
    "Mortality, 2nd+ infection Juveniles, Winter",
    "Mortality, 2nd+ infection Adults, Winter",
    "Mortality, 2nd+ infection Adults, Summer",
    "Mortality, Recovered Juveniles, Summer",
    "Mortality, Recovered Adults, Summer",
    "Mortality, Recovered Juveniles, Winter",
    "Mortality, Recovered Adults, Winter",
    "Beta, Susceptible Juveniles, Winter",
    "Beta, Susceptible Juveniles, Summer",
    "Beta, Recovered Adults, Winter",
    "Beta, Recovered Juveniles, Winter",
    "Beta, Recovered Adults, Summer",
    "Beta, Recovered Juveniles, Summer",
    "Recovery, 1st infection Juveniles, Summer",
    "Recovery, 1st infection Adults, Summer",
    "Recovery, 2nd+ infection Juveniles, Summer",
    "Rcevory, 2nd+ infection Adults, Summer",
    "Recovery, 1st infection Juveniles, Winter",
    "Recovery, 1st infection Adults, Winter",
    "Recovery, 2nd+ infection Juveniles, Winter",
    "Recovery, 2nd+ infection Adults, Winter"
  )
  # Remove values in each column that are outside the range of values in
  # the same column in the prior samples
  posterior_samples <- abc$adj.values |> as.data.frame() |>
    mutate(Distribution = "Posterior") |>
    mutate(across(all_of(param_names), ~ ifelse(. < min(prior_samples[[cur_column()]]),
                                                 min(prior_samples[[cur_column()]]),
                                                 ifelse(. > max(prior_samples[[cur_column()]]),
                                                        max(prior_samples[[cur_column()]]), .))))


  # Create a data frame for plotting
  plot_data <- bind_rows(posterior_samples, prior_samples) |>
    select(-c("mortality_Sa_summer", "mortality_Ra_summer")) |>
    pivot_longer(cols = !any_of("Distribution"), names_to = "Parameter",
                 values_to = "Value") |>
    mutate(Parameter = factor(Parameter, levels = param_names,
           labels = param_labels))

  # Create faceted plot
  tidyplot(plot_data, x = Value, fill = Distribution) |>
    add(geom_density(alpha = 0.7)) |>
    adjust_x_axis_title("Parameter Value") |>
    adjust_title("Prior vs Posterior Distributions") |>
    split_plot(Parameter, ncol = 2, nrow = 2)
}
