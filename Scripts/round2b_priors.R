library(data.table)
library(here)
library(abc)
library(triangle)
library(qs2)
library(mc2d)
library(fitdistrplus)

options(scipen = 999)

demo_params <- fread(here("Data/Input/sample_data_round1.csv"))
abc_first_pass <- qs_read(here("Data/Validation/abc_round1_2variable.qs2"))
round1_metrics <- fread(here("Data/Validation/round1_validation_metrics.csv"))

demo_params_selected <- demo_params[round1_metrics$dc, ][, `:=`(
  Region = abc_first_pass$region,
  Weights = 1 /
    (abc_first_pass$dist +
      .Machine$double.eps)
)][Region == TRUE, ]
summary(demo_params_selected)

demo_params_selected[1:10, ]
