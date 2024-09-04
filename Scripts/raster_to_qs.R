library(purrr)
library(stringr)
library(qs)
library(terra)
library(doParallel)
library(foreach)

# Set up
num_workers <- 4
output_dir <- "/glade/work/pilowskyj/Round1_matrix"
input_dir <- "/glade/work/pilowskyj/scratch"
incomplete_sims <- read.delim("/glade/work/pilowskyj/Round1/incomplete_simulations.txt", sep = "\n")
cl <- makeCluster(num_workers)
registerDoParallel(cl)


input_subdir <- system("cd /glade/work/pilowskyj/scratch; ls -f", intern = TRUE) |>
  _[-c(1:2)] |> discard(\(string) str_detect(string, "txt"))

input_files <- map(input_subdir, function(dir) {
  path <- file.path(input_dir, dir)
  system(paste0("cd ", path, "; ls -f"), intern = TRUE) |> _[-c(1:2)]
})

process_raster <- function(path, i, j, lookup = NULL) {
  if (is.null(lookup)) {
    output_subdir <- file.path(output_dir, input_subdir[i])
  } else {
    index <- input_subdir[i] |> str_extract("[0-9]+") |> as.integer()
    output_subdir <- file.path(output_dir,
                               paste0("simulation", incomplete_sims[index,]))
  }
  output_path <- file.path(output_subdir,
                           paste0(input_files[[i]][j] |>
                                    str_split_i("\\.", 1), ".qs"))
  input_path <- file.path(input_dir, input_subdir[i], path)
  if (!dir.exists(output_subdir)) {
    dir.create(output_subdir)
  }
  if (!file.exists(output_path)) {
    rast(input_path) |> as.array() |>
      qsave(output_path)
  }
}

# Parallel processing of the rasters
foreach(i = seq_along(input_files),
        .packages = c("purrr", "stringr", "qs", "terra")) %dopar% {
  iwalk(input_files[[i]], function(path, j) {
    process_raster(path, i, j, lookup = incomplete_sims)
  })
}

# Stop the cluster after processing is complete
stopCluster(cl)

