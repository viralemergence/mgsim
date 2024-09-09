library(tidyverse)
library(epizootic)
set.seed(198)
runs5000 <- replicate(5000, aspatial_siri(
  initial_pop = c(50000, 50000, 0, 1, 0, 0, 0, 0),
  season_length = 100,
  mortality = c(0.004, 0, 0.00505, 0.00105, 0.004, 0, 0.0045, 5e-04),
  fecundity = 15/182,
  transmission = c(0.00002, 0.00001, 7.84e-06, 3.92e-06),
  recovery = c(0.05714286, 0.05714286, 0.1, 0.1),
  carrying_capacity = 150000,
  abundance_threshold = 10,
  season = "breeding"
))
set.seed(198)
seasonal5000 <- replicate(5000, aspatial_siri_seasonal(
  initial_pop = c(50000, 50000, 0, 1, 0, 0, 0, 0),
  mortality = c(0.004, 0, 0.00505, 0.00105, 0.004, 0, 0.0045, 5e-04) |>
    (`*`)(100),
  fecundity = 0.08241758 |> (`*`)(100),
  transmission = c(0.00002, 0.00001, 7.84e-06, 3.92e-06) |>
    (`*`)(100),
  recovery = rep(1, 4),
  carrying_capacity = 150000,
  abundance_threshold = 10,
  season = "breeding"
))
list(runs5000, seasonal5000) |> map(t) |> map(as.data.frame) |>
  map(\(x) set_names(x, c("Sj", "Sa", "I1j", "I1a", "Rj", "Ra",
                                                        "I2j", "I2a"))) |>
  bind_rows() |>
  mutate(Simulator = c(rep("Daily", 5000), rep("Seasonal", 5000))) |>
  pivot_longer(Sj:I2a, names_to = "State", values_to = "Count") |>
  ggplot(aes(y = Count, x = State, fill = State)) + geom_boxplot() +
  scale_fill_paletteer_d(`"dichromat::Categorical_12"`, direction = -1) +
  facet_wrap(~Simulator) +
  scale_y_continuous(trans = "log1p")
