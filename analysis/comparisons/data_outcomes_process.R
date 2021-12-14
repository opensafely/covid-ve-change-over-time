################################################################################

# This script:


################################################################################

library(tidyverse)
library(glue)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobs <- FALSE
  group <- "02"
  
} else{
  removeobs <- TRUE
  group <- args[[1]]
}

study_parameters <- readr::read_rds(here::here("output", "lib", "study_parameters.rds"))

# outcome data after the start of comparison 1 
# (for unvax with 2 brands, use min start date)
# as this is the earliest possible date that we will look for outcome events
data_in <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds")) %>%
  filter(comparison == 1) %>%
  distinct(patient_id, time_zero_date, end_fu_date) %>%
  mutate(end_fu_date = end_fu_date + (study_parameters$n_comparisons-1)*28) %>%
  # keep only the earliest date for each individual
  arrange(patient_id, time_zero_date) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  select(patient_id, 
         earliest_date = time_zero_date,
         latest_date = end_fu_date)

################################################################################
## join outcomes data
data_outcomes <- data_in %>%
  left_join(
    readr::read_rds(here::here("output", "data", "data_long_postest_dates.rds")) %>%
      select(patient_id, postest_date = date),
    by = "patient_id"
  ) %>%
  left_join(
    readr::read_rds(here::here("output", "data", "data_long_covidadmitted_dates.rds")) %>%
      select(patient_id, covidadmitted_date = date),
    by = "patient_id"
  ) %>%
  left_join(
    arrow::read_feather(
      here::here("output", "input_covs.feather")) %>%
      select(patient_id, coviddeath_date, death_date, dereg_date) %>%
      mutate(across(ends_with("date"), as.Date)),
    by = "patient_id"
  ) %>%
  # in case coviddeath_date and death_date different dates
  mutate(across(c(coviddeath_date, death_date),
                ~ if_else(
                  !is.na(coviddeath_date) & !is.na(death_date),
                  min(coviddeath_date, death_date),
                  .x
                ))) %>%
  # add outcome for noncoviddeath
  mutate(
    noncoviddeath_date = if_else(
      !is.na(death_date) & is.na(coviddeath_date),
      death_date, 
      as.Date(NA_character_))
  ) %>%
  mutate(across(c(postest_date, covidadmitted_date, coviddeath_date, death_date, noncoviddeath_date),
                ~ if_else(!is.na(.x) & earliest_date < .x & .x <= latest_date,
                          .x,
                          as.Date(NA_character_)))) %>%
  select(-earliest_date, -latest_date)
  

readr::write_rds(
  data_outcomes,
  here::here("output", glue("jcvi_group_{group}"), "data", "data_outcomes.rds"),
  compress = "gz")