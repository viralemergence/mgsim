library(here)
library(tidyverse)
library(qs)
library(poems)
region <- qread(here("Data/Input/finch_region.qs"))
daily_sims <- here("Data/Output/Round 1 daily") %>%
  list.files(full.names = T) %>%
  map(qread) %>%
  map("abundance")
seasonal_sims <- here("Data/Output/Round 1 daily") %>%
  list.files() %>%
  str_extract("[0-9]+") %>%
  map(\(x) qread(here("Data/Output/Round 1 seasonal",
                      paste0("sample_", x, "_results.qs")))) %>%
  map("abundance")
map2(daily_sims, seasonal_sims, t.test, paired = TRUE) # Not all differences
# are in the same direction. Sometimes the daily sim population is higher,
# sometimes the reverse. Observe!
region$raster_from_values(daily_sims[[9]][, 59, 2] - seasonal_sims[[9]][ , 59, 2]) %>%
  raster::plot(main = "Delta between daily and seasonal simulation abundance")
region$raster_from_values(daily_sims[[10]][, 59, 2] - seasonal_sims[[10]][ , 59, 2]) %>%
  raster::plot(main = "Delta between daily and seasonal simulation abundance")
