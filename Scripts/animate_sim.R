animate_sim <- function(path, region, burn_in, years) {

  sim <- qread(path)
  arr <- sim$abundance[, burn_in:ncol(sim$abundance),]
  season <- c(0, 0.5)
  # Create a data.frame with the correct x, y, z coordinates
  df <- expand.grid(Population = 1:dim(arr)[1], y = 1:dim(arr)[2],
                    z = 1:dim(arr)[3])

  # Convert the data.frame to a data.table
  dt <- as.data.table(df)

  # Add the corresponding values from the array
  dt[, Abundance := arr[cbind(Population, y, z)]]
  dt$Year <- years[dt$y]
  dt$Season <- season[dt$z]
  dt$x <- region$coordinates[df$Population, 1]
  dt$y <- region$coordinates[df$Population, 2]
  dt$Time <- dt$Year + dt$Season

  # Plot the data
  anim <- ggplot() +
    geom_sf(data = basemap) +
    geom_tile(data = dt, mapping = aes(x, y, fill = Abundance)) +
    scale_fill_viridis_c(labels = label_comma()) +
    transition_time(Time) +
    labs(title = "Year: {frame_time}") +
    theme_void()

  animate(anim, nframes = length(unique(dt$Time)))
}
