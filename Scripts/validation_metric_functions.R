library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(qs)
library(Rcpp)
#### Validation metrics functions ####

##### Fill missing years in a data frame ####
fill_missing_years <- function(df, start_year, end_year) {
  # Use NA or a special marker to indicate missing/zero years instead of using a file path
  missing_marker <- NA_character_

  # Full sequence of years
  all_years <- seq(start_year, end_year)

  # Identify missing years
  missing_years <- setdiff(all_years, df$Year)

  # If there are no missing years, return the dataframe as is
  if (length(missing_years) == 0) {
    return(df)
  }

  # Identify a template year that is present in the data
  template_year <- df[df$Year == min(df$Year), ]  # Assuming the first year as the template

  # Function to create rows for a missing year based on the template
  create_missing_year_rows <- function(year) {
    template <- template_year
    template$Year <- year
    template$path <- missing_marker  # Mark the missing years with NA or another marker
    return(template)
  }

  # Create dataframes for each missing year
  missing_dfs <- lapply(missing_years, create_missing_year_rows)

  # Combine the original dataframe with the missing years dataframes
  df_filled <- bind_rows(c(list(df), missing_dfs))

  # Order by Year (and any other columns if needed)
  df_filled <- df_filled[order(df_filled$Year), ]

  return(df_filled)
}

##### Convert flat coordinate to a 2d coordinate in a matrix #####
convert_flat_to_2d <- function(flat_index, nrows) {
  row <- (flat_index - 1) %% nrows + 1
  col <- (flat_index - 1) %/% nrows + 1
  return(c(row, col))
}

###### Function to either read the file or return the zero array #####

# Define an array of zeroes to use for "crashed" simulations
zero_array <- array(0, dim = c(106, 161, 1))

read_or_zero <- function(path) {
  if (is.na(path)) {  # Check if path is NA (missing year)
    return(zero_array)
  } else {
    return(qread(path))  # Otherwise, read the file
  }
}

##### Helper function to sum abundances at BCR indices for house finch trends #####
cppFunction('
NumericVector sum_positions(NumericVector sim_array, IntegerMatrix positions, int n_years, int n_rows, int n_cols) {
    // Create a vector to store the total abundance for each year
    NumericVector total_abundance(n_years);

    // Traverse the positions matrix and sum over the relevant parts of the array
    for (int k = 0; k < n_years; k++) {
        double sum = 0.0;
        for (int i = 0; i < positions.nrow(); i++) {
            int row = positions(i, 0) - 1;  // Adjust for zero-indexing in C++
            int col = positions(i, 1) - 1;
            sum += sim_array(row + col * n_rows + k * n_rows * n_cols);  // Index into the 3D array
        }
        total_abundance[k] = sum;
    }
    return total_abundance;
}
')

##### Penalty function for house finch trends #####
abundance_trend_penalty <- function(percent_change, trend_upper, trend_lower) {
  if (any(is.na(c(trend_upper, trend_lower)))) {
    stop("The trend inputs must not be NA")
  }
  if (is.na(percent_change)) {
    penalty <- 20
  } else if (percent_change < trend_lower) {
    penalty <- trend_lower - percent_change
  } else if (percent_change > trend_upper) {
    penalty <- percent_change - trend_upper
  } else {
    penalty <- 0
  }
  return(penalty)
}

##### Calculate abundance trend metric #####

# Optimized trend metrics calculation using Rcpp
calculate_trend_metrics <- function(year_range, baseline_year) {
  n_years <- length(1970:2016)
  results <- data.frame(
    Abundance = numeric(n_years * 33),
    bcr = integer(n_years * 33),
    Year = integer(n_years * 33)
  )

  # Track index for filling preallocated data frame
  row_index <- 1

  # Efficient subsetting and summing using C++ with Rcpp
  for (int_val in 5:37) {
    positions <- which(bcr == int_val, arr.ind = TRUE)

    if (nrow(positions) > 0) {
      # Call the C++ function for summation
      total_abundance <- sum_positions(flat_sim_array, positions, n_years, 106, 161)

      # Fill the preallocated result
      results[row_index:(row_index + n_years - 1), ] <- data.frame(
        Abundance = total_abundance,
        bcr = rep(int_val, n_years),
        Year = 1970:2016
      )

      row_index <- row_index + n_years  # Update index for next BCR
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

  baseline <- map_dbl(trend_by_bcr$model, \(m) predict(m, newdata = data.frame(Year = baseline_year)))
  trend_by_bcr$percent_change <- (trend_by_bcr$estimate / baseline) * 100

  penalty <- trend_by_bcr |>
    select(BCR = bcr, percent_change) |>
    right_join(if (baseline_year == 1970) trend1970 else trend1993, by = "BCR") |>
    mutate(penalty = abundance_trend_penalty(percent_change, estimate_ucl, estimate_lcl)) |>
    pull(penalty) |>
    sum()

  return(penalty)
}

##### Penalty function for house finch arrival dates #####
hm_arrival_function <- function(observed, low, high) {

  if (is.na(observed)) {
    return(76)
  }

  if (observed < low) {
    penalty <- low - observed
  } else if (observed > high) {
    penalty <- observed - high
  } else {
    penalty <- 0
  }

  return(penalty)
}

##### Penalty function for Mycoplasma arrival dates #####
# Penalty function
mg_arrival_function <- function(observed, low, high) {

  if (is.na(observed)) {
    return(54)
  }

  if (observed < low) {
    penalty <- low - observed
  } else if (observed > high) {
    penalty <- observed - high
  } else {
    penalty <- 0
  }

  return(penalty)
}

##### Penalty function for spatiotemporal point prevalence #####
prevalence_penalty <- function(prevalence, lower, upper, zero) {
  if (is.na(prevalence)) {
    return(0.25)
  }

  if (!is.na(zero)) {
    penalty <- plogis(prevalence) - plogis(0)
  } else if (prevalence < lower) {
    penalty <- plogis(lower) - plogis(prevalence)
  } else if (prevalence > upper) {
    penalty <- plogis(prevalence) - plogis(upper)
  } else {
    penalty <- 0
  }

  return(penalty)
}
