# process tests data
library(tidyverse)
library(lubridate)

data_tests <- arrow::read_feather(
  file = here::here("output", "input_tests.feather")) %>%
  mutate(across(c(contains("_date")), 
                ~ floor_date(
                  as.Date(., format="%Y-%m-%d"),
                  unit = "days"))) %>%
  select(-elig_date)

readr::write_rds(
  data_tests,
  here::here("output", "data", "data_tests.rds"),
  compress = "gz"
)
