library(tidyverse)
library(lubridate)
library(glue)

dummy_data <- arrow::read_feather(here::here("analysis", "lib", "dummy_data_vax.feather"))


dummy_data_covs <- dummy_data %>%
  select(patient_id, age_1, age_2, sex, elig_date, region_0) %>%
  mutate(
    # fill in these vars one study_definition_covs reviewed
  )

arrow::write_feather(dummy_data_covs, here::here("analysis", "covs", "dummy_data_covs.feather"))
  