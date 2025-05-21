library(stringr)
library(here)
library(readr)
library(doParallel)
library(foreach)
i_am("mgsim/Scripts/validation_debug.R")
source("/glade/u/home/pilowskyj/mgsim/Scripts/validation_metric_functions.R")
# Register parallel backend with doParallel
num_cores <- 120
cl <- makeCluster(num_cores, type = "FORK", outfile = "")
registerDoParallel(cl)
results_dir <- system(
  "cd /glade/work/pilowskyj/Round1_matrix; ls -f",
  intern = T
) |>
  gtools::mixedsort() |>
  _[-c(1:2)]
round1_priors <- read_csv(here(
  "mgsim/Data_minimal/Input/sample_data_round1.csv"
))
year_lookup <- data.frame(index = 1:77, Year = 1940:2016)
dc <- read_csv(
  "/glade/u/home/pilowskyj/mgsim/Data_minimal/Validation/dc.csv"
) |>
  _$dc
# Load the presence/absence data
presabs <- here(
  "mgsim/Data_minimal/Validation/haemorhous_presence_absence.csv"
) |>
  read_csv()

# Gather relevant data
presence_list <- results_dir[dc] |>
  map(\(f) paste0("cd /glade/work/pilowskyj/Round1_matrix/", f)) |>
  lapply(list.files, full.names = TRUE) |>
  lapply(function(p) {
    data.frame(
      path = p,
      index = as.numeric(str_split_i(p, "/", 8) |> str_split_i("_", 1)),
      infected = str_detect(p, "I")
    ) |>
      arrange(index) |>
      left_join(year_lookup, by = join_by(index))
  })

# Check presence_list
print(paste("Length of presence_list:", length(presence_list)))
flush.console()

# Check the first element of presence_list (if not too large)
if (length(presence_list) > 0) {
  print("First element of presence_list:")
  print(head(presence_list[[1]]))
  flush.console()
}

presabs_metric <- foreach(
  plist = presence_list,
  .combine = c,
  .packages = c("dplyr", "purrr", "qs")
) %dopar%
  {
    tryCatch(
      {
        # Print the start of the loop
        print(paste("Starting loop for presence list element at:", Sys.time()))
        flush.console()

        # Original code

        plist_filled <- plist |>
          fill_missing_years(1940, 2016) |>
          group_split(Year)
        sim_array <- array(0, dim = c(106, 161, length(plist_filled)))

        # Adding together the matrices for all age classes and disease compartments for each year
        for (i in seq_along(plist_filled)) {
          df <- plist_filled[[i]]
          combined_files <- Reduce(`+`, lapply(df$path, read_or_zero))
          sim_array[,, i] <- combined_files
        }

        # Calculating presence/absence penalty
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
          ) |>
          pull(penalty) |>
          sum()

        # Print successful completion of loop
        print(paste("Completed loop successfully at:", Sys.time()))
        flush.console()

        return(penalty)
      },
      error = function(e) {
        # Print error message
        print(paste("Error in loop at:", Sys.time(), "-", conditionMessage(e)))
        flush.console()
        return(NA)
      }
    )
  }

presabs <- data.frame(sim = c(1:10000)[dc], hm_presabs = presabs_metric)
write_csv(
  presabs,
  "/glade/u/home/pilowskyj/mgsim/Data_minimal/Validation/presabs.csv"
)
