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
  right_join(data_eligible_e %>% select(patient_id, start_1_date, min_elig_date), 
            by = "patient_id") %>%
  mutate(across(patient_id, as.integer)) %>%
  mutate(across(ends_with("_date"), as.Date))

nrows <- nrow(dummy_data)

dates_seq <- seq(as.Date("2021-01-01"), as.Date("2021-11-30"), 1)

K <- study_parameters$max_comparisons

# function for result_test_k_date
test_k_date <- function(k, test_result = "any") {
  
  name <- glue("{test_result}_test_{k}_date")
  
  date <- sample(dates_seq, size = nrows, replace = TRUE)
  
  dummy_data %>%
    mutate(!! sym(name) := date) %>%
    mutate(across(!! sym(name),
                  ~ case_when(
                    start_1_date + days((k-1)*28) <  .x & .x <= start_1_date + days(k*28) ~ .x,
                    TRUE ~ as.Date(NA_character_)
                  ))) %>%
    select(!! sym(name))
  
}


# function for result_test_k_date so that any_test_k_n >= positive_test_k_n
test_k_n <- function(k) {
  
  name_any <- glue("any_test_{k}_n")
  name_positive <- glue("positive_test_{k}_n")
  
  dummy_data %>%
    mutate(!! sym(name_positive) := rpois(n = nrow(.), lambda = 0.25)) %>%
    mutate(!! sym(name_any) := !! sym(name_positive) + rpois(n = nrow(.), lambda = 1)) %>%
    select(!! sym(name_any), !! sym(name_positive))
  
}

# pregnancy dates
preg_k_date <- function(k, type = "pregnancy") {
  
  if (type == "pregnancy") {
    name <- glue("preg_36wks_{k}_date")
  } else {
    name <- glue("pregdel_pre_{k}_date")
  }
  
  date <- sample(dates_seq, size = nrows, replace = TRUE)
  
  dummy_data %>%
    transmute(!! sym(name) := date) 
  
}
# pregnancy indicators
preg_k <- function(k) {
  
  name <- glue("preg_{k}")
  
  dummy_data %>%
    transmute(!! sym(name) := as.integer(rbernoulli(n = nrow(.), p=0.01)))
  
}


dummy_data_tests <- dummy_data %>%
  bind_cols(lapply(1:(K+1), test_k_date)) %>%
  bind_cols(lapply(1:(K+1), test_k_n)) %>%
  bind_cols(lapply(1:K, preg_k_date)) %>%
  bind_cols(lapply(1:K, function(x) preg_k_date(k=x, type = "delivery"))) %>%
  bind_cols(lapply(1:K, preg_k)) %>%
  mutate(test_hist_1_n = rpois(n = nrow(.), lambda = 3),
         test_hist_2_n = rpois(n = nrow(.), lambda = 3),
         test_hist_3_n = rpois(n = nrow(.), lambda = 3)) %>%
  mutate(across(ends_with("date"), as.POSIXct))


arrow::write_feather(dummy_data_tests, here::here("analysis", "dummy_data_tests.feather"))
  