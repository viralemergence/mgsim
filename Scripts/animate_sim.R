library(gganimate)
library(sp)
library(sf)
library(data.table)

animate_sim <- function(array, region, years = 1940:2016, burn_in = 0) {
  if (burn_in > 0) {
    arr <- array[,, burn_in:dim(array)[3]]
    timesteps <- dim(arr)[3]
  } else {
    arr <- array
    timesteps <- 154
  }
  coords <- region$region_raster |> coordinates()
  # Create a data.frame with the correct x, y, z coordinates
  df <- expand.grid(y = 1:dim(arr)[2], x = 1:dim(arr)[1], z = 1:dim(arr)[3])

  # Convert the data.frame to a data.table
  dt <- as.data.table(df)

  # Add the corresponding values from the array
  dt[, Abundance := arr[cbind(x, y, z)]]
  dt$Year <- years[round(dt$z / 2) + 1]
  dt$Season <- if_else(dt$z %% 2 == 0, 0, 0.5)
  dt$x <- rep(coords[, 1], timesteps)
  dt$y <- rep(coords[, 2], timesteps)
  dt$Time <- dt$Year + dt$Season
  dt$Abundance <- if_else(dt$Abundance == 0, NA_real_, dt$Abundance)

  # Create a basemap
  basemap <- rnaturalearth::ne_coastline() |>
    st_transform(st_crs(region$region_raster)) |>
    st_crop(st_bbox(region$region_raster))

  # Plot the data
  anim <- ggplot() +
    geom_sf(data = basemap) +
    geom_tile(data = dt, mapping = aes(x, y, fill = Abundance)) +
    scale_fill_viridis_c(
      labels = scales::label_comma(),
      na.value = "transparent"
    ) +
    transition_time(Time) +
    labs(title = "Year: {round(frame_time)}") +
    theme_void()

  return(gganimate::animate(
    anim,
    nframes = length(unique(dt$Time)),
    fps = 4,
    render = gifski_renderer()
  ))
}
