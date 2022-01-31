######################################

# This script:
# - creates dummy data for study_definition_tests.py

######################################

library(tidyverse)
library(lubridate)
library(glue)

# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

source(here::here("analysis", "lib", "dummy_data_functions.R"))

# individuals eligible based on box d criteria
data_eligible_e <- readr::read_csv(
  here::here("output", "data", "data_eligible_e.csv"))

set.seed(5476)

dummy_data <- arrow::read_feather(
  file = here::here("analysis", "dummy_data.feather")) %>%
  select(patient_id, elig_date) %>%
  right_join(data_eligible_e %>% select(patient_id, start_1_date), 
            by = "patient_id") 

nrows <- nrow(dummy_data)

dates_seq <- seq(as.Date("2021-01-01"), as.Date("2021-11-30"), 1)

K <- study_parameters$max_comparisons

# function for any_test_k_date
test_k_n <- function(k, test_result = "any") {
  
  name <- glue("{test_result}_test_{k}_date")
  
  date <- sample(dates_seq, size = nrows, replace = TRUE)
  
  dummy_data %>%
    mutate(!! sym(name) := as.POSIXct(date)) %>%
    mutate(across(!! sym(name),
                  ~ case_when(
                    start_1_date + days((k-1)*28) <  .x & .x <= start_1_date + days(k*28) ~ .x,
                    TRUE ~ NA_POSIXct_
                  ))) %>%
    select(!! sym(name))
  
}



dummy_data %>%
  bind_cols(lapply(1:(K+1), test_k_n))



  