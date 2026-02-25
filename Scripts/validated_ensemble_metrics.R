library(ggplot2)
library(stringr)
library(here)
library(readr)
library(furrr)
plan(multicore)
source(here("Scripts/validation_metric_functions.R"))

results_dir <- here("Data/Output/Round3a") |>
  list.files(full.names = T) |>
  gtools::mixedsort()
year_lookup <- data.frame(index = 1:77, Year = 1940:2016)
filter_paths <- function(file_paths) {
  # Use str_extract to capture the last number before the underscore
  last_numbers <- as.numeric(str_extract(file_paths, "(?<=/)(\\d+)(?=_)"))

  # Filter paths where the last number is between 63 and 77
  filtered_paths <- file_paths[last_numbers >= 63 & last_numbers <= 77]

  return(filtered_paths)
}

# Put file paths in the right order
file_paths <- results_dir |>
  map(list.files, full.names = T) |>
  map(filter_paths)
data_list <- map(file_paths, \(fp) {
  data.frame(
    path = fp,
    season = str_extract(fp, "summer|winter"),
    number = as.integer(str_extract(fp, "(?<=/)(\\d+)(?=_)")),
    infected = str_detect(fp, "I")
  ) |>
    mutate(season = factor(season, levels = c("winter", "summer"))) |>
    arrange(number, season)
})
simulation_id <- results_dir |> str_sub(67) |> as.integer()
process_single_number <- function(df, number) {
  df %>%
    filter(number == !!number) %>%
    group_split(infected) %>%
    map(\(subset_df) pull(subset_df, path)) %>%
    map(\(l) map(l, qread)) |>
    map(\(l) map(l, \(m) m[mg_trends$region_cell])) |>
    map(\(df) Reduce(`+`, df)) |>
    setNames(c("healthy", "sick")) |>
    bind_cols() |>
    mutate(index = number, region_cell = mg_trends$region_cell)
}

# Load the presence/absence data
pres_hfds <- read_csv(here(
  "Data_minimal/Validation/mycoplasma_presence.csv"
))

# Extract the file paths for infected birds
infected_list <- results_dir |>
  map(list.files, full.names = TRUE) |>
  lapply(function(p) {
    data.frame(
      path = p,
      index = as.numeric(str_split_i(p, "/", 10) |> str_split_i("_", 1)),
      infected = str_detect(p, "I")
    ) |>
      filter(infected) |>
      arrange(index) |>
      left_join(year_lookup, by = join_by(index)) |>
      filter(Year %in% 1994:2016)
  })

no_infection <- sapply(infected_list, function(d) nrow(d) == 0)

mg_presence_df <- future_map_dfr(
  infected_list[!no_infection],
  \(df) {
    # Fill missing years and group by Year
    df_filled <- fill_missing_years(df, 1994, 2016) |> group_split(Year)

    # Pre-allocate a 3D array for the simulation (106x161xyears)
    sim_array <- array(0, dim = c(106, 161, length(df_filled)))

    # Efficiently sum matrices into the 3D array, one slice per year
    for (i in seq_along(df_filled)) {
      year_df <- df_filled[[i]]

      # Load and sum the matrices for all files in this year
      combined_files <- Reduce(`+`, lapply(year_df$path, read_or_zero))

      # Assign the combined matrix to the appropriate slice
      sim_array[,, i] <- combined_files
    }

    # Calculate the metric for this simulation
    metric <- pres_hfds |>
      rowwise() |>
      mutate(
        rowcol = list(convert_flat_to_2d(region_cell, 106)),
        sick = sum(sim_array[rowcol[[1]], rowcol[[2]], min_index:max_index]),
        none_sick = sick == 0
      )
    return(metric)
  },
  .id = "simulation_id"
)

# Load the presence/absence data
presabs <- here(
  "Data_minimal/Validation/haemorhous_presence_absence.csv"
) |>
  read_csv()

# Gather relevant data
presence_list <- results_dir |>
  lapply(list.files, full.names = TRUE) |>
  lapply(function(p) {
    data.frame(
      path = p,
      index = as.numeric(str_split_i(p, "/", 10) |> str_split_i("_", 1)),
      infected = str_detect(p, "I")
    ) |>
      arrange(index) |>
      left_join(year_lookup, by = join_by(index))
  })

hm_presabs_df <- future_map_dfr(
  presence_list,
  \(plist) {
    # Filling in years where the simulation did not write because the population was zero with paths to an array of all zeroes
    plist_filled <- plist |> fill_missing_years(1940, 2016) |> group_split(Year)

    # This 3D array represents a yearly simulation from 1940 to 2016, with time as the z axis
    sim_array <- array(0, dim = c(106, 161, length(plist_filled)))

    # Adding together the matrices for all age classes and disease compartments for each year
    for (i in seq_along(plist_filled)) {
      df <- plist_filled[[i]]

      combined_files <- Reduce(`+`, lapply(df$path, read_or_zero))

      sim_array[,, i] <- combined_files
    }

    penalty <- presabs |>
      mutate(
        rowcol = map(region_cell, ~ convert_flat_to_2d(.x, 106)),
        occupancy = pmap_dbl(
          list(rowcol, min_index, max_index),
          function(rc, min_i, max_i) {
            sum(sim_array[rc[1], rc[2], min_i:max_i] > 0)
          }
        ),
        penalty = pmap_dbl(
          list(always_present, occupancy, min_index, max_index),
          function(always_present_val, occ, min_i, max_i) {
            if (always_present_val) {
              length(min_i:max_i) - occ
            } else {
              occ
            }
          }
        )
      )
    return(penalty)
  },
  .id = "simulation_id"
)

# Read in trend data
trend1993 <- read_csv(here(
  "Data_minimal/Validation/abundance_trend_1993on.csv"
))
trend1970 <- read_csv(here(
  "Data_minimal/Validation/abundance_trend_1970on.csv"
))
# and conservation regions
bcr <- here("Data_minimal/Validation/bird_conservation_regions.qs") |>
  qread()

# Gather relevant data
abundance_list <- results_dir |>
  lapply(list.files, full.names = TRUE) |>
  lapply(function(p) {
    data.frame(
      path = p,
      index = as.numeric(str_split_i(p, "/", 10) |> str_split_i("_", 1)),
      infected = str_detect(p, "I")
    ) |>
      arrange(index) |>
      left_join(year_lookup, by = join_by(index)) |>
      filter(Year %in% 1970:2016)
  })

trend_df <- future_map_dfr(
  abundance_list,
  \(alist) {
    # Preprocess data and fill missing years
    alist_filled <- alist |> fill_missing_years(1970, 2016) |> group_split(Year)

    # Preallocate a 3D array for sim_array results
    sim_array_results <- array(0, dim = c(106, 161, length(alist_filled)))

    # Process each year chunk
    for (i in seq_along(alist_filled)) {
      df <- alist_filled[[i]]

      # Read and combine files
      paths <- df$path
      files <- lapply(paths, read_or_zero)
      combined_files <- Reduce(`+`, files)

      # Fill the array
      sim_array_results[,, i] <- combined_files
    }

    # Flatten sim_array_results for use in C++
    flat_sim_array <- as.vector(sim_array_results)

    year_range <- 1993:2016
    baseline_year <- 1970
    n_years <- length(1970:2016)
    results <- data.frame(
      Abundance = numeric(n_years * 33),
      bcr = integer(n_years * 33),
      Year = integer(n_years * 33)
    )

    # Track index for filling preallocated data frame
    row_index <- 1

    for (int_val in 5:37) {
      positions <- which(bcr == int_val, arr.ind = TRUE)

      if (nrow(positions) > 0) {
        total_abundance <- sum_positions(
          flat_sim_array,
          positions,
          n_years,
          106,
          161
        )

        # Fill the preallocated result
        results[row_index:(row_index + n_years - 1), ] <- data.frame(
          Abundance = total_abundance,
          bcr = rep(int_val, n_years),
          Year = 1970:2016
        )

        row_index <- row_index + n_years # Update index for next BCR
      }
    }

    results <- results[1:(row_index - 1), ]

    # Continue with linear modeling and penalty calculation
    trend_by_bcr <- results |>
      filter(Year %in% year_range) |>
      group_by(bcr) |>
      nest() |>
      mutate(model = map(data, \(df) lm(Abundance ~ Year, data = df))) |>
      mutate(tidy = map(model, tidy)) |>
      unnest(tidy) |>
      filter(term == "Year")

    baseline <- map_dbl(
      trend_by_bcr$model,
      \(m) predict(m, newdata = data.frame(Year = baseline_year))
    )
    trend_by_bcr$percent_change <- (trend_by_bcr$estimate / baseline) * 100
    penalty1970 <- trend_by_bcr |>
      select(BCR = bcr, percent_change) |>
      right_join(
        if (baseline_year == 1970) trend1970 else trend1993,
        by = "BCR"
      ) |>
      mutate(
        penalty = abundance_trend_penalty(
          percent_change,
          estimate_ucl,
          estimate_lcl
        ),
        trend = baseline_year
      )

    baseline_year <- 1993
    year_range <- 1993:2016

    trend_by_bcr <- results |>
      filter(Year %in% year_range) |>
      group_by(bcr) |>
      nest() |>
      mutate(model = map(data, \(df) lm(Abundance ~ Year, data = df))) |>
      mutate(tidy = map(model, tidy)) |>
      unnest(tidy) |>
      filter(term == "Year")

    baseline <- map_dbl(
      trend_by_bcr$model,
      \(m) predict(m, newdata = data.frame(Year = baseline_year))
    )
    trend_by_bcr$percent_change <- (trend_by_bcr$estimate / baseline) * 100
    penalty1993 <- trend_by_bcr |>
      select(BCR = bcr, percent_change) |>
      right_join(
        if (baseline_year == 1970) trend1970 else trend1993,
        by = "BCR"
      ) |>
      mutate(
        penalty = abundance_trend_penalty(
          percent_change,
          estimate_ucl,
          estimate_lcl
        ),
        trend = baseline_year
      )
    penalty <- bind_rows(penalty1970, penalty1993)
    return(penalty)
  },
  .id = "simulation_id"
)

winter_indices <- seq(1, 39, 2)
summer_indices <- seq(2, 40, 2)
# Read in validation data
prevalence <- here(
  "Data_minimal/Validation/mycoplasma_point_prevalence.csv"
) |>
  read_csv() |>
  left_join(year_lookup) |>
  rowwise() |>
  mutate(
    index = if_else(
      Season == "Summer",
      summer_indices[index - 55],
      winter_indices[index - 55]
    )
  )

prevalence_list <- results_dir |>
  lapply(list.files, full.names = TRUE) |>
  lapply(function(p) {
    data.frame(
      path = p,
      index = as.numeric(str_split_i(p, "/", 10) |> str_split_i("_", 1)),
      season = str_extract(p, "summer|winter"),
      infected = str_detect(p, "I")
    ) |>
      mutate(season = factor(season, levels = c("winter", "summer"))) |>
      arrange(index, season) |>
      left_join(year_lookup, by = join_by(index)) |>
      filter(Year %in% 1995:2014)
  })

point_prevalence_metric <- future_map_dfr(
  prevalence_list,
  \(plist) {
    plist_filled <- plist |>
      fill_missing_years(start_year = 1995, end_year = 2014) |>
      group_split(Year, season)

    # Preallocate a 3D array for sim_array results
    sim_array_results <- array(0, dim = c(106, 161, length(plist_filled)))

    # Process each year chunk to manage memory and avoid redundant computations
    for (i in seq_along(plist_filled)) {
      df <- plist_filled[[i]]

      # Read files only once per chunk
      paths <- df$path
      infected <- lapply(paths[1:4], read_or_zero) |> Reduce(f = `+`, x = _)
      total <- lapply(paths, read_or_zero) |> Reduce(f = `+`, x = _)
      prev <- infected / total

      # Fill the preallocated array
      sim_array_results[,, i] <- prev
    }

    # Calculate the metric for this presence_list
    penalty <- prevalence |>
      mutate(
        rowcol = list(convert_flat_to_2d(region_cell, 106)),
        prevalence = sim_array_results[rowcol[1], rowcol[2], index],
        penalty = prevalence_penalty(prevalence, Lower, Upper, Zero)
      )

    gc()
    return(penalty)
  },
  .id = "simulation_id"
)

mg_arrival <- here(
  "Data_minimal/Validation/mycoplasma_first_arrival.csv"
) |>
  read_csv() |>
  group_by(region_cell) |>
  summarize(
    first_arrival = first(first_arrival) |> round(),
    first_arrival_high = first(first_arrival_high) |> round(),
    first_arrival_low = first(first_arrival_low) |> round()
  ) |>
  left_join(year_lookup, by = c("first_arrival_high" = "Year")) |>
  left_join(
    year_lookup,
    by = c("first_arrival_low" = "Year"),
    suffix = c("_high", "_low")
  )

mg_arrival_metric <- future_map_dfr(
  infected_list[!no_infection],
  \(df) {
    # Fill missing years and group by Year
    df_filled <- fill_missing_years(df, 1994, 2016) |> group_split(Year)

    # Pre-allocate a 3D array for the simulation (106x161xyears)
    sim_array <- array(0, dim = c(106, 161, length(df_filled)))

    # Efficiently sum matrices into the 3D array, one slice per year
    for (i in seq_along(df_filled)) {
      year_df <- df_filled[[i]]

      # Sum the matrices for all files in this year
      combined_files <- Reduce(`+`, lapply(year_df$path, read_or_zero))

      # Assign the summed matrix to the appropriate slice
      sim_array[,, i] <- combined_files
    }

    arrival_index <- map(
      mg_arrival$region_cell,
      ~ convert_flat_to_2d(.x, 106)
    ) |>
      map_int(function(rc) {
        presences <- which(sim_array[rc[1], rc[2], ] > 0)
        if (length(presences) == 0) {
          NA_integer_
        } else {
          min(presences, na.rm = TRUE)
        }
      })

    penalty <- mg_arrival |>
      mutate(arrival_index = arrival_index) |>
      rowwise() |>
      mutate(
        penalty = mg_arrival_function(arrival_index, index_low, index_high)
      )

    gc()
    return(penalty)
  },
  .id = "simulation_id"
)

presence_list <- results_dir |>
  lapply(list.files, full.names = TRUE) |>
  lapply(function(p) {
    data.frame(
      path = p,
      index = as.numeric(str_split_i(p, "/", 10) |> str_split_i("_", 1)),
      infected = str_detect(p, "I")
    ) |>
      arrange(index) |>
      left_join(year_lookup, by = join_by(index))
  })

hm_arrival <- here(
  "Data_minimal/Validation/haemorhous_first_arrival.csv"
) |>
  read_csv() |>
  filter(region_cell != 6075) |>
  group_by(region_cell) |>
  summarize(
    first_arrival = first(first_arrival) |> round(),
    first_arrival_high = first(first_arrival_high) |> round(),
    first_arrival_low = first(first_arrival_low) |> round()
  ) |>
  left_join(year_lookup, by = c("first_arrival_high" = "Year")) |>
  left_join(
    year_lookup,
    by = c("first_arrival_low" = "Year"),
    suffix = c("_high", "_low")
  )

hm_arrival_df <- future_map_dfr(
  presence_list,
  \(df) {
    # Fill missing years and group by Year
    df_filled <- fill_missing_years(df, 1940, 2016) |> group_split(Year)

    # Pre-allocate a 3D array for the simulation (106x161xyears)
    sim_array <- array(0, dim = c(106, 161, length(df_filled)))

    # Efficiently sum matrices into the 3D array, one slice per year
    for (i in seq_along(df_filled)) {
      year_df <- df_filled[[i]]

      # Sum the matrices for all files in this year
      combined_files <- Reduce(`+`, lapply(year_df$path, read_or_zero))

      # Assign the summed matrix to the appropriate slice
      sim_array[,, i] <- combined_files
    }

    # Compute the arrival index for each region
    arrival_index <- map(
      hm_arrival$region_cell,
      ~ convert_flat_to_2d(.x, 106)
    ) |>
      map_int(function(rc) {
        presences <- which(sim_array[rc[1], rc[2], ] > 0)
        if (length(presences) == 0) {
          NA_integer_
        } else {
          min(presences, na.rm = TRUE)
        }
      })

    # Calculate penalties based on the arrival index
    penalty <- hm_arrival |>
      mutate(arrival_index = arrival_index) |>
      rowwise() |>
      mutate(
        penalty = hm_arrival_function(arrival_index, index_low, index_high)
      )

    gc()
    return(penalty)
  },
  .id = "simulation_id"
)
