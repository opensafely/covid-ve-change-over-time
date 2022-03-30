################################################################################

# This script:
# - calculates min and max follow-up dates for each subgroup

################################################################################

library(tidyverse)

# create output directories ----
fs::dir_create(here::here("output", "lib"))

# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

## read data
# covariates data
data_covariates <- readr::read_rds(
  here::here("output", "data", "data_covariates.rds")) 

# processed data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

data_min_max_fu <- data_covariates %>%
  inner_join(
    data_processed,
    by = "patient_id"
  ) %>%
  group_by(subgroup) %>%
  summarise(
    min_fu_date = min(start_k_date),
    max_fu_date = max(end_k_date),
    .groups = "keep"
  ) %>% 
  ungroup() %>%
  mutate(across(max_fu_date,
                ~ pmin(as.Date(study_parameters$end_date), .x)))

readr::write_rds(
  data_min_max_fu,
  here::here("output", "lib", "data_min_max_fu.rds")
)
