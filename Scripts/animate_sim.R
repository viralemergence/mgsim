library(gganimate)
library(biscale)
library(sp)
library(sf)
library(data.table)

animate_sim <- function(
  array,
  region,
  years = 1940:2016,
  burn_in = 0,
  remove_outliers = FALSE
) {
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

  # Handle outliers if requested
  if (remove_outliers) {
    threshold_95 <- quantile(dt$Abundance, 0.95, na.rm = TRUE)
    dt$Abundance <- pmin(dt$Abundance, threshold_95)
  }

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
      na.value = "transparent",
      oob = scales::squish
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

prevalence_sim <- function(sim, region, data_dir, biclass = FALSE) {
  # Assertions
  assert_integerish(
    sim,
    lower = 1,
    upper = 10000,
    any.missing = FALSE,
    len = 1
  )
  assert_class(region, classes = "Region")
  assert_directory_exists(data_dir, access = "r")

  # Load the data
  dir <- file.path(data_dir, paste0("simulation", sim))
  assert_directory_exists(dir, access = "r")
  num_slashes <- str_count(dir, "/")
  year_lookup <- data.frame(index = 1:77, Year = 1940:2016)
  file_info <- dir |>
    list.files(full.names = TRUE) |>
    lapply(
      function(p) {
        data.frame(
          path = p,
          index = as.numeric(
            str_split_i(p, "/", num_slashes + 2) |>
              str_split_i("_", 1)
          ),
          season = str_extract(p, "summer|winter")
        ) |>
          mutate(
            season = factor(
              season,
              levels = c(
                "winter",
                "summer"
              )
            )
          ) |>
          arrange(index, season) |>
          left_join(year_lookup, by = join_by(index)) |>
          filter(str_detect(path, "Sj|Sa|I1j|I1a|Rj|Ra|I2j|I2a"))
      }
    ) |>
    bind_rows() |>
    fill_missing_years(start_year = 1940, end_year = 2016) |>
    filter(Year > 1992)

  coords <- region$region_raster |> coordinates()
  prevalence_df <- map(seq_len(nrow(file_info)), function(i) {
    # Detect the compartment and life stage from the file name
    compartment <- str_extract(file_info$path[i], "Sj|Sa|I1j|I1a|Rj|Ra|I2j|I2a")
    life_stage <- if_else(str_detect(compartment, "j"), "Juvenile", "Adult")
    compartment <- str_remove(compartment, "j|a")
    # Read the .qs file
    arr <- read_or_zero(file_info$path[i])
    # Create a data.frame with the correct x, y, z coordinates
    df <- expand.grid(y = 1:dim(arr)[2], x = 1:dim(arr)[1])

    # Convert the data.frame to a data.table
    dt <- as.data.table(df)

    # Add the corresponding values from the array
    dt[, Abundance := arr[cbind(x, y)]]
    dt$Year <- file_info$Year[i]
    dt$Season <- file_info$season[i]
    dt$x <- coords[, 1]
    dt$y <- coords[, 2]
    dt$Time <- if_else(dt$Season == "winter", dt$Year, dt$Year + 0.5)
    dt$Compartment <- compartment
    dt$Life_Stage <- life_stage

    return(dt)
  }) |>
    bind_rows() |>
    mutate(Infected = str_detect(Compartment, "I")) |>
    group_by(Infected, Time, x, y) |>
    summarize(Abundance = sum(Abundance, na.rm = T), .groups = "drop") |>
    pivot_wider(
      names_from = Infected,
      values_from = Abundance,
      names_prefix = "Infected"
    ) |>
    mutate(
      Prevalence = InfectedTRUE / (InfectedTRUE + InfectedFALSE),
      Population_Size = InfectedTRUE + InfectedFALSE
    )

  if (biclass) {
    prevalence_df <- bi_class(
      prevalence_df,
      Prevalence,
      Population_Size,
      style = "fisher",
      dim = 4
    )
  }

  # Create a basemap
  basemap <- rnaturalearth::ne_coastline() |>
    st_transform(st_crs(region$region_raster)) |>
    st_crop(st_bbox(region$region_raster))

  # Create the legend once
  if (biclass) {
    legend <- bi_legend(
      pal = "GrPink2",
      dim = 4,
      xlab = "Prevalence",
      ylab = "Population Size",
      size = 8
    )
  }

  if (biclass) {
    p <- ggplot() +
      geom_sf(data = basemap) +
      geom_tile(
        data = prevalence_df,
        mapping = aes(x, y, fill = bi_class),
        show.legend = FALSE
      ) +
      bi_scale_fill(
        pal = "GrPink2",
        dim = 4,
        na.value = "transparent"
      ) +
      annotation_custom(
        grob = ggplotGrob(legend),
        xmin = 2000000,
        xmax = 3700000,
        ymin = -1300000,
        ymax = 400000
      ) +
      transition_time(Time) +
      labs(title = "Year: {round(frame_time)}") +
      theme_void()
  } else {
    p <- ggplot() +
      geom_sf(data = basemap) +
      geom_tile(
        data = prevalence_df,
        mapping = aes(x, y, fill = Prevalence),
        show.legend = TRUE
      ) +
      scale_fill_viridis_c(
        na.value = "transparent",
        oob = scales::squish
      ) +
      transition_time(Time) +
      labs(title = "Year: {round(frame_time)}") +
      theme_void()
  }

  return(gganimate::animate(
    p,
    nframes = length(unique(prevalence_df$Time)),
    fps = 4,
    render = gifski_renderer()
  ))
}
