library(qs)
library(tidyverse)
library(epizootic)
library(furrr)
plan(multisession)
data_dir <- here::here("Data/Output/Round 1")
region <- here::here("Data/Input/finch_region.qs") |> qread()
dc <- list.files(data_dir, pattern = "qs", full.names = TRUE) |>
  gtools::mixedsort() %>%
  future_map_lgl(function(file) {
    mat <- qread(file)
    mat$abundance[3531, 59, 2] > 0
  })
ranges <- list.files(data_dir, pattern = "qs", full.names = TRUE) |>
  gtools::mixedsort() %>%
  future_map_int(function(file) {
    mat <- qread(file)
    sum(mat$abundance[, 59, 2] > 0)
  })
sim <- list.files(data_dir, pattern = "qs", full.names = TRUE) |>
  gtools::mixedsort() %>% .[[44]] |> qread()
region$raster_from_values(sim$abundance[ , 59, 2]) |> raster::plot(xlim = c(0, 3e06))
sim$abundance[3531, ,]
