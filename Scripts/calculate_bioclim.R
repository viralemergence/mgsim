library(terra)
library(tidyterra)
library(fs)
library(tidyverse)
library(furrr)
proj_crs <- "+proj=aea +lon_0=-94.5 +lat_1=21.5 +lat_2=47.5 +lat_0=34.5 +datum=WGS84 +units=m +no_defs"
urban_rast <- here::here("Data/Input/urban_raster_latlong.grd")
directory <- "/Volumes/MEDIA/CHELSA"
predict_string <- "~/Documents/Very_Large_Data/predictors/CHELSA/predictors_"
finch_region <- qs::qread(here::here("Data/Input/finch_region.qs"))
process_data <- function(y) {
  library(terra)
  urban <- urban_rast |> rast()
  gap_filler <- here::here("Data/Input/region_mask.grd") |> rast()
  climlist <- directory |> dir_ls() |>
    map(dir_ls, regexp = y) |>
    map(gtools::mixedsort) |>
    map(rast) |>
    map(subst, from = -32768, to = NA) |>
    map(aggregate, fact = 10) |>
    map(crop, y = urban) |>
    map(project, gap_filler) |>
    map(mask, gap_filler) |>
    set_names(nm = c("prec", "tmax", "tmin")) |>
    identity()
  urban_rast_proj <- urban |>
    focal(fun = mean, na.rm = TRUE, weights = matrix(1, 3, 3)) |>
    project(gap_filler) |>
    mask(gap_filler)
  bios <- bioclima::clima(c(1:19), tmin = climlist$tmin,
                          tmax = climlist$tmax,
                          prcp = climlist$prec) |>
    c(urban_rast_proj[[paste0("AD", y)]])
  writeRaster(bios, paste0(predict_string, y, ".grd"))
}
plan(multisession, workers = 3)
c(1940:2016) %>%
  discard(function(y) {
    paste0(predict_string, y, ".grd") |> file.exists()
  }) %>%
  future_walk(.f = process_data,
              .progress = T)

c(1940:1949) %>%
  future_walk(function(y) {
    gap_filler <- here::here("Data/Input/region_mask.grd") |> rast()
    urban <- urban_rast |> rast() %>% .[[paste0("AD", y)]] |>
      focal(fun = mean, na.rm = TRUE, weights = matrix(1, 3, 3)) |>
      project(gap_filler) |>
      mask(gap_filler)
    fix <- paste0(predict_string, y, ".grd") |> rast() %>%
      .[[1:19]] |> c(urban)
    r <- fix + 0
    writeRaster(r, paste0(predict_string, y, ".grd"), overwrite = T)
  }, .progress = T)
