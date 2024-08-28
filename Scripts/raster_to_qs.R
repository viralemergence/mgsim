library(tidyverse)
library(qs)
library(here)
library(furrr)
plan(multisession, workers = 4)
output_dir <- here("Data/Output/Round1_matrix")
input_dir <- here("Data/Output/Round1")
input_subdir <- list.dirs(input_dir) |>
  gtools::mixedsort() |> _[-1] |>
  keep(\(l) length(list.files(l)) > 0)
input_files <- map(input_subdir, list.files, full.names = T)
process_raster <- function(path, i, j) {
  output_subdir <- file.path(output_dir, str_split_i(input_subdir[i], "/", 9))
  output_path <- file.path(output_subdir,
                           paste0(str_split_i(input_files[[i]][j], "/", 10) |>
                                    str_split_i("\\.", 1), ".qs"))
  if (!dir.exists(output_subdir)) {
    dir.create(output_subdir)
  }
  if (!file.exists(output_path)) {
    rast(path) |> as.array() |>
      qsave(output_path)
  }
}
input_files |>
  future_iwalk(\(paths, i) iwalk(paths, \(path, j) process_raster(path, i, j)))
