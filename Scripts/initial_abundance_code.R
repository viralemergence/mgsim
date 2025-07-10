library(terra)
library(poems)
library(qs)
# Create a mask for the native range of the house finch
data_dir <- "/Users/caryinstitute/Documents/mgsim/Data/Input"
template <- rast("Data/Input/habitat_suitability_v1.tif")
native_range <- rast("~/Downloads/native range2_modified.tif")
native_range[native_range < 180] <- NA
native_range[native_range >= 180] <- 1
native_range_mask <- native_range |> _[[4]] |> extend(template) |>
  resample(template, "sum")
native_range_mask[native_range_mask > 1] <- 1
native_range_mask[native_range_mask < 1] <- NA
native_range_mask |> writeRaster("Data/Input/native_range_mask.tif")
# Geocode Jones Beach
c(40.595833, -73.515278) |> rev() |> matrix(nrow = 1) |>
  project(crs("EPSG:4326"), crs(template)) |> cellFromXY()
region <- data_dir |> file.path("finch_region.qs") |> qread()
which(region$region_indices == 6075)
