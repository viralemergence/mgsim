library(qs)
library(tidyverse)
library(gganimate)
library(paletteer)
library(cowplot)
library(here)
library(abind)
library(data.table)
library(sf)
sim <- 12

# Load data
round1 <- paste0("Data/Output/Round 1 daily/sample_", sim, "_results.qs") %>%
  here() |> qread()
round2 <- paste0("Data/Output/Round 2 daily/sample_", sim, "_results.qs") %>%
  here() |> qread()
region <- here("Data/Input/finch_region.qs") |> qread()
years <- c(1940:2016)
basemap <- rnaturalearth::ne_coastline() |>
  st_transform(st_crs(region$region_raster)) |>
  st_crop(st_bbox(region$region_raster))

#### Abundance visualization ####
# Assemble abundance
full_sim_abundance <- abind(round1$abundance[, 6:58, ],
                            round2$abundance[, 1:23, ], along = 2)
season <- c(0, 0.5)
# Create a data.frame with the correct x, y, z coordinates
df <- expand.grid(Population = 1:dim(full_sim_abundance)[1],
                  y = 1:dim(full_sim_abundance)[2],
                  z = 1:dim(full_sim_abundance)[3])

# Convert the data.frame to a data.table
dt <- as.data.table(df)

# Add the corresponding values from the array
dt[, Abundance := full_sim_abundance[cbind(Population, y, z)]]
dt$Year <- years[dt$y]
dt$Season <- season[dt$z]
dt$x <- region$coordinates[df$Population, 1]
dt$y <- region$coordinates[df$Population, 2]
dt$Time <- dt$Year + dt$Season

# Plot the data
anim <- ggplot() +
  geom_sf(data = basemap) +
  geom_tile(data = dt, mapping = aes(x, y, fill = Abundance)) +
  scale_fill_viridis_c(labels = scales::label_comma()) +
  transition_time(Time) +
  labs(title = "Year: {frame_time}") +
  theme_void()

gganimate::animate(anim, nframes = length(unique(dt$Time)))

#### Presence + prevalence visualization ####
infected <- round2$abundance_segments$stage_1_compartment_2 +
  round2$abundance_segments$stage_2_compartment_2 +
  round2$abundance_segments$stage_1_compartment_4 +
  round2$abundance_segments$stage_2_compartment_4
filler <- array(0, dim = c(6355, 53, 2))
full_sim_infected <- abind(filler, infected, along = 2)
df <- expand.grid(Infected = 1:dim(full_sim_infected)[1],
                  y = 1:dim(full_sim_infected)[2],
                  z = 1:dim(full_sim_infected)[3])
dt[, Presence := Abundance > 0]
dt[, Infected := full_sim_infected[cbind(df$Infected, df$y, df$z)]]
dt[, Prevalence := Infected / Abundance]

anim2 <- ggplot() +
  geom_sf(data = basemap) +
  geom_tile(data = dt, mapping = aes(x, y, fill = Prevalence)) +
  scale_fill_viridis_c(labels = scales::label_comma()) +
  transition_time(Time) +
  labs(title = "Year: {frame_time}") +
  theme_void()

gganimate::animate(anim2, nframes = length(unique(dt$Time)))
