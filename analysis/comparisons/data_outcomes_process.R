######################################

# What this script does:
# imports data created by the `data_comparisons_process.R` script
# converts hospitalisation and infection episodes into long format
# saves as a one-row-per-event dataset

######################################

# Preliminaries ----

## Import libraries ----
library(tidyverse)
library(glue)

# Import custom user functions from lib
# source(here("lib", "utility_functions.R"))

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

## Import data_comparisons ----
data_outcomes <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_covariates.rds")) %>%
  distinct(patient_id) %>%
  left_join(
    arrow::read_feather(here::here("output", "input_covs.feather")) %>%
      select(patient_id, 
             matches(c("^positive\\_test\\_\\d+\\_date",
                       "^primary\\_care\\_covid_case\\_\\d+\\_date",
                       "^primary\\_care\\_suspected\\_covid\\_\\d+\\_date",
                       "^covidadmitted_\\d+\\_date",
                       "^admitted\\_unplanned\\_infectious\\_\\d+\\_date",
                       "^discharged\\_unplanned\\_infectious\\_\\d+\\_date",
                       "^admitted\\_unplanned\\_\\d+\\_date",
                       "^discharged\\_unplanned\\_\\d+\\_date")),
             coviddeath_date, death_date),
             by = "patient_id") %>%
  mutate(across(ends_with("date"), as.Date))


## create one-row-per-event datasets ----
# for positive test, hospitalisation/discharge, covid in primary care, death

data_admissions <- data_outcomes %>%
  select(patient_id, 
         matches("^admitted\\_unplanned\\_\\d+\\_date"),
         matches("^discharged\\_unplanned\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(".value", "index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_drop_na = TRUE
  ) %>%
  select(patient_id, 
         index, 
         admitted_date=admitted_unplanned, 
         discharged_date = discharged_unplanned) %>%
  arrange(patient_id, admitted_date)

data_admissions_infectious <- data_outcomes %>%
  select(patient_id, 
         matches("^admitted\\_unplanned\\_infectious\\_\\d+\\_date"), 
         matches("^discharged\\_unplanned\\_infectious\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(".value", "index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_drop_na = TRUE
  ) %>%
  select(patient_id,
         index,
         admitted_date=admitted_unplanned_infectious, 
         discharged_date = discharged_unplanned_infectious) %>%
  arrange(patient_id, admitted_date)

#remove infectious admissions from all admissions data
data_admissions_noninfectious <- anti_join(
  data_admissions,
  data_admissions_infectious,
  by = c("patient_id", "admitted_date", "discharged_date")
)


data_pr_suspected_covid <- data_outcomes %>%
  select(patient_id,
         matches("^primary\\_care\\_suspected\\_covid\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "suspected_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

data_pr_probable_covid <- data_outcomes %>%
  select(patient_id,
         matches("^primary\\_care\\_covid\\_case\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "probable_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

data_postest <- data_outcomes %>%
  select(patient_id, 
         matches("^positive\\_test\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "postest_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)



readr::write_rds(
  data_admissions, 
  here::here("output", group,  "data", "data_long_admission_dates.rds"), 
  compress="gz")
readr::write_rds(
  data_admissions_infectious, 
  here::here("output", group, "data", "data_long_admission_infectious_dates.rds"), 
  compress="gz")
readr::write_rds(
  data_admissions_noninfectious, 
  here::here("output", group, "data", "data_long_admission_noninfectious_dates.rds"), 
  compress="gz")
readr::write_rds(
  data_pr_probable_covid, 
  here::here("output", group, "data", "data_long_pr_probable_covid_dates.rds"),
  compress="gz")
readr::write_rds(
  data_pr_suspected_covid, 
  here::here("output", group, "data", "data_long_pr_suspected_covid_dates.rds"), 
  compress="gz")
readr::write_rds(
  data_postest, 
  here::here("output", group,  "data", "data_long_postest_dates.rds"), 
  compress="gz")
