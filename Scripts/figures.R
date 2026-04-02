#### Setup ####
# Load packages
library(conflicted)
library(tidyverse)
library(tidyplots)
library(tidymodels)
library(DALEXtra)
library(qs)
library(here)
library(paletteer)
library(furrr)
library(scales)
library(cowplot)
library(magick)
library(rasterVis)
library(terra)
conflicts_prefer(dplyr::filter, dplyr::select, dplyr::lag)
# Source other scripts as needed
source(here("Scripts/validation_metric_functions.R"))
source(here("Scripts/convenience_functions.R"))
source(here("Scripts/animate_sim.R"))
# Geospatial data
region <- qread(here("Data/Input/finch_region.qs"))
# Validation metrics
round2_metrics <- read_csv(here(
  "Data/Validation/round2c_validation_metrics.csv"
))
trend_metrics2 <- read_csv(here("Data/Validation/round2c_trend_metrics.csv")) |>
  pivot_longer(
    everything(),
    names_to = c("metric", "sim"),
    names_sep = "\\.\\.\\."
  )
round2_metrics$hm_trend_1970[round2_metrics$dc] <- trend_metrics2 |>
  filter(metric == "penalty1970") |>
  pull(value)
round2_metrics$hm_trend_1993[round2_metrics$dc] <- trend_metrics2 |>
  filter(metric == "penalty1993") |>
  pull(value)
round1_metrics <- read_csv(here(
  "Data/Validation/round1_validation_metrics.csv"
)) |>
  mutate(hm_trend_1970 = rep(NA_real_), hm_trend_1993 = rep(NA_real_))
round1_metrics$hm_trend_1970[round1_metrics$dc] <- round1_metrics$hm_trend[seq(
  1,
  9752,
  2
)]
round1_metrics$hm_trend_1993[round1_metrics$dc] <- round1_metrics$hm_trend[seq(
  2,
  9752,
  2
)]
round3a_metrics <- read_csv(here(
  "Data/Validation/round3a_validation_metrics.csv"
))
trend_metrics3a <- read_csv(here(
  "Data/Validation/round3a_trend_metrics.csv"
)) |>
  pivot_longer(
    everything(),
    names_to = c("metric", "sim"),
    names_sep = "\\.\\.\\."
  )
round3a_metrics$hm_trend_1970[round3a_metrics$dc] <- trend_metrics3a |>
  filter(metric == "penalty1970") |>
  pull(value)
round3a_metrics$hm_trend_1993[round3a_metrics$dc] <- trend_metrics3a |>
  filter(metric == "penalty1993") |>
  pull(value)
# Priors
round1_priors <- read_csv(here("Data/Input/sample_data_round1.csv"))
round3a_priors <- read_csv(
  here("Data/Input/sample_data_round3a.csv"),
  show_col_types = FALSE
) |>
  mutate(
    mortality_Sa_summer = 0,
    mortality_Ra_summer = mortality_Sa_summer
  )
# Selected models
abc31 <- qread(here("Data/Validation/abc_round3a.qs"))
selected_samples_3 <- round3a_priors |>
  rowid_to_column("sample") |>
  semi_join(abc31$unadj.values, copy = TRUE) |>
  pull(sample)
selected_metrics <- round3a_metrics[selected_samples_3, ]

#### Figure 1: Validation Metrics ####
validation_metrics <- bind_rows(
  round1_metrics,
  round2_metrics,
  round3a_metrics,
  selected_metrics
) |>
  mutate(
    Round = factor(c(
      rep("1", 10000),
      rep("2", 10000),
      rep("3", 10000),
      rep("Selected", 32)
    ))
  ) |>
  select(-hm_trend, -index)
p1 <- validation_metrics |>
  filter(Round != "Selected") |>
  tidyplot(x = Round, fill = dc) |>
  add_barstack_absolute() |>
  adjust_y_axis_title("Count") |>
  adjust_x_axis_title("Simulation Round") |>
  adjust_colors(new_colors = paletteer_d("nbapalettes::nets")) |>
  adjust_legend_title("Finches in DC in 1994")
p2 <- validation_metrics |>
  pivot_longer(
    c(mg_presence, hm_arrival),
    values_to = "Value",
    names_to = "Metric"
  ) |>
  filter(!is.na(Value)) |>
  mutate(
    Metric = recode(
      Metric,
      mg_presence = "MG persistence",
      hm_arrival = "Finch arrival"
    )
  ) |>
  tidyplot(x = Value, fill = Round) |>
  add_histogram(bins = 30, alpha = 0.5, position = "identity") |>
  adjust_y_axis_title("Count") |>
  adjust_colors(new_colors = paletteer_d("PNWColors::Sunset2", 4)) |>
  split_plot(by = Metric, scales = "free_x") |>
  adjust_y_axis(
    transform = transform_pseudo_log(),
    force_continuous = TRUE,
    breaks = c(0, 10, 100, 1000),
    labels = label_number(big.mark = ",", accuracy = 1)
  )

tmp_p1 <- tempfile(fileext = ".png")
tmp_p2 <- tempfile(fileext = ".png")

ggplot2::ggsave(
  filename = tmp_p1,
  plot = p1,
  width = 8,
  height = 2.4,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = tmp_p2,
  plot = p2,
  width = 8,
  height = 4.6,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

img1 <- magick::image_read(tmp_p1) |> magick::image_trim()
img2 <- magick::image_read(tmp_p2) |> magick::image_trim()
spacer <- magick::image_blank(
  width = max(magick::image_info(img1)$width, magick::image_info(img2)$width),
  height = 40,
  color = "white"
)

fig1 <- c(
  img1,
  spacer,
  img2
) |>
  magick::image_append(stack = TRUE)

magick::image_write(
  fig1,
  path = here("Visualizations/figure1_POM.png"),
  format = "png"
)

#### Figure 2: Selected Posterior Distributions ####
# Placeholder. This will be a plot of selected posterior distributions compared
# to the priors.

#### Figure 3: First arrival dates of MG ####

results <- process_sims(
  selected_samples_3,
  abc31$weights,
  here("Data/Output/Round3a")
)
prevalence <- results$prevalence

# Check the 3D prevalence array to see which time slice is the first to show a
# prevalence > 0 at each XY location
first_arrival <- apply(prevalence, c(1, 2), function(x) {
  if (all(x == 0, na.rm = TRUE)) {
    return(NA) # No arrival
  } else {
    years <- seq(1940, 2016.5, by = 0.5) |> round()
    index <- which(x > 0)[1]
    return(years[index]) # First time slice with prevalence > 0
  }
})
first_arrival_raster <- region$raster_from_values(t(first_arrival)[
  region$region_indices
])
# Calculate the weighted centroid of the prevalence distribution for each time slice
weightedCentre <- function(x, y, z) {
  require(matrixStats) #; require(Hmisc)
  if (anyNA(c(x, y, z))) {
    stop("There are missing values present in x, y, or z")
  }
  if (length(z) != length(x)) {
    stop("Number of weights supplied not equal to number of coordinates")
  }
  x_y = cbind(x, y)
  w = z
  mean = matrixStats::colWeightedMeans(x_y, w = w, na.rm = FALSE)
  wght_cent <- matrix(mean, ncol = 2) |> as.data.frame()
  colnames(wght_cent) <- c("x", "y")
  return(wght_cent)
}

centroids <- map_dfr(seq_len(dim(prevalence)[3]), function(slice) {
  x <- region$coordinates$x
  y <- region$coordinates$y
  z <- as.vector(t(prevalence[,, slice])[region$region_indices])
  z[is.na(z)] <- 0 # Replace NA values with 0 for weighting
  weightedCentre(x, y, z)
}) |>
  mutate(year = seq(1940, 2016.5, by = 0.5)) |>
  mutate(distance_west_lag1 = lag(x, n = 1) - x)

fig3 <- gplot(first_arrival_raster) +
  geom_tile(aes(fill = value)) +
  scale_fill_paletteer_c(
    palette = "pals::parula",
    name = "First Arrival Year",
    na.value = "grey70"
  ) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "right") +
  geom_point(
    data = centroids,
    aes(x = x, y = y),
    color = "black",
    size = 2,
    shape = 21,
    fill = "white"
  )
ggplot2::ggsave(
  filename = here("Visualizations/figure3_mg_arrival.png"),
  plot = fig3,
  width = 6,
  height = 6,
  dpi = 300
)

#### Figure 4: Regression model results ####
# Prep data
abc31_processed <- process_sims(
  selected_samples_3,
  abc31$weights,
  here("Data/Output/Round3a")
)
total_population_size <- abc31_processed$total_population
total_cases <- prevalence * total_population_size
arr <- total_cases[,, 110:154]
arr2 <- total_population_size[,, 110:154]
years <- 1994:2016
coords <- region$region_raster |> raster::coordinates()
dt <- expand.grid(y = 1:dim(arr)[2], x = 1:dim(arr)[1], z = 1:dim(arr)[3]) |>
  mutate(Cases = arr[cbind(x, y, z)], Population_size = arr2[cbind(x, y, z)])
dt$Year <- years[round(dt$z / 2) + 1]
dt$Season <- if_else(dt$z %% 2 == 0, 0, 0.5)
dt$x <- rep(coords[, 1], 45)
dt$y <- rep(coords[, 2], 45)
dt$Time <- dt$Year + dt$Season
index_case_coords <- c(y = 621173.8, x = 1475944)
dt$Distance_from_index_case <- sqrt(
  (dt$x - index_case_coords["x"])^2 + (dt$y - index_case_coords["y"])^2
)

# Extract environmental predictors
col_names <- 1:19 %>%
  as.character() %>%
  map(~ paste0("bio", .)) %>%
  append("urban") %>%
  flatten_chr()
files2 <- list.files(
  "~/Documents/Very_Large_Data/predictors/CHELSA",
  pattern = "^predictors_.+",
  full.names = T
) %>%
  keep(str_detect, pattern = "\\.grd$") %>%
  .[55:77] |>
  set_names(years)
extracted <- dt |>
  group_by(Year) |>
  group_split() |>
  map(~ select(.x, x, y)) |>
  future_map2(
    files2,
    function(p, f) {
      r <- rast(f)
      terra::extract(r, p, ID = F)
    },
    .options = furrr_options(seed = T)
  ) %>%
  map(~ set_names(., col_names)) %>%
  list_rbind(names_to = "Year") |>
  select(bio1, bio4, urban) |>
  bind_cols(dt)

# Fit model
analysis_input <- extracted |>
  select(
    bio1,
    bio4,
    urban,
    Population_size,
    Cases,
    Year,
    Distance_from_index_case,
    Season
  ) |>
  mutate(Season = if_else(Season %% 0.5 == 0, "Winter", "Summer")) |>
  mutate(Season = factor(Season, levels = c("Winter", "Summer"))) |>
  drop_na()
urb_split <- initial_split(analysis_input, strata = "Cases")
urb_train <- training(urb_split)
urb_test <- testing(urb_split)
urb_recipe <- recipe(Cases ~ ., data = urb_train) |>
  step_normalize(-Season, -all_outcomes()) |>
  step_dummy(Season)
tune_spec <- boost_tree(
  # tuned already to save time
) |>
  set_engine("xgboost") |>
  set_mode("regression")
# urb_folds <- vfold_cv(urb_train, v = 5, strata = "Cases")
# urb_tune <- workflow() |>
#     add_recipe(urb_recipe) |>
#     add_model(tune_spec) |>
#     tune_grid(
#         resamples = urb_folds,
#         grid = 20,
#         control = control_grid(save_pred = TRUE)
#     )
best_tree <- list(
  trees = 316,
  tree_depth = 4,
  min_n = 12,
  loss_reduction = 0.120,
  sample_size = 0.621,
  mtry = 8,
  learn_rate = 0.210
)
final_wf_urb <- workflow() |>
  add_model(tune_spec) |>
  add_recipe(urb_recipe) |>
  finalize_workflow(best_tree)
final_fit_urb <- final_wf_urb |> last_fit(urb_split)

# Plot
urb_explainer <- explain_tidymodels(
  final_fit_urb$.workflow[[1]],
  y = urb_test$Cases,
  data = urb_test |>
    select(
      bio1,
      bio4,
      urban,
      Year,
      Distance_from_index_case,
      Season,
      Population_size
    ),
  label = ""
)
label_map <- c(
  bio1 = "Mean annual temperature (bio1)",
  bio4 = "Temperature seasonality (bio4)",
  urban = "Urban land use",
  Population_size = "Population size",
  Year = "Year",
  Distance_from_index_case = "Distance from index case"
)

var_imp <- urb_explainer |> variable_importance()
urb_profile <- urb_explainer %>%
  model_profile(variables = "urban", N = 2000)
fig4a <- plot(var_imp) +
  labs(subtitle = NULL, caption = NULL) +
  scale_x_discrete(labels = \(x) dplyr::recode(x, !!!label_map, .default = x))
fig4b <- plot(urb_profile) +
  labs(subtitle = NULL, caption = NULL) +
  scale_y_continuous(name = "Average predicted number of cases") +
  scale_x_continuous(name = "Urbanization index") +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Partial dependence profile: Urbanization")

fig4 <- plot_grid(
  fig4a,
  fig4b,
  ncol = 1,
  rel_heights = c(1.2, 1)
)

ggplot2::ggsave(
  filename = here("Visualizations/figure4_regression.png"),
  plot = fig4,
  width = 6,
  height = 8,
  dpi = 300
)

#### Figure S1 ####
s1 <- plot_abc_posteriors(abc31, round1_priors)
tidyplots::save_plot(s1, here("Visualizations/Figure_S1_posteriors.pdf"))
