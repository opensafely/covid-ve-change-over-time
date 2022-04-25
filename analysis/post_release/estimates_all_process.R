library(tidyverse)

release_folder <- "release_20220401"

# read data
estimates_all <- readr::read_csv(
  here::here(release_folder, "estimates_all.csv"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))

# output folder
fs::dir_create(here::here("output", "metareg", "data"))

# prcoess
data_metareg <- estimates_all %>%
  filter(
    !reference_row,
    variable %in% "k",
    model %in% "unadjusted2" # error in the labeling, this corresponds to unadjusted
    ) %>%
  mutate(model = "adjusted") %>%
  select(subgroup, comparison, outcome, model, k = label, estimate, conf.low, conf.high) %>%
  mutate(across(subgroup, factor, levels = 1:4, labels = subgroups))

# save data
readr::write_csv(
  data_metareg,
  here::here(release_folder, "data_metareg.csv"))
