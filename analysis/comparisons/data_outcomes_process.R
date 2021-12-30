################################################################################

# This script:
# read ids and earliest start_fu_date from all eligible individuals
# create data_outcomes: one-row-per-individual of outcome data occurring after the
# earliest time zero
# save data_outcomes
################################################################################

library(tidyverse)
library(glue)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  group <- "02"
  
} else{
  group <- args[[1]]
}

# read comparisons data
data_ids <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds")) %>%
  # only keep patient_id and earliest start_fu_date
  # Don't just filter to comparison 1, as unvax arm has two rows per comparison,
  # one for each brand. Just keep the ealiest.
  arrange(patient_id, start_fu_date) %>%
  select(patient_id, start_fu_date) %>%
  distinct(patient_id, .keep_all = TRUE)

################################################################################
## join outcomes data
data_outcomes <- data_ids %>%
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
    readr::read_rds(
      here::here("output", "data", "data_covs.rds")) %>%
      select(patient_id, coviddeath_date, death_date, dereg_date) %>%
      mutate(across(ends_with("date"), as.Date)),
    by = "patient_id"
  ) %>%
  # in case coviddeath_date and death_date different dates
  mutate(across(c(coviddeath_date, death_date),
                ~ if_else(
                  !is.na(coviddeath_date) & !is.na(death_date),
                  pmin(coviddeath_date, death_date, na.rm = TRUE),
                  .x
                ))) %>%
  # add outcome for noncoviddeath
  mutate(
    noncoviddeath_date = if_else(
      !is.na(death_date) & is.na(coviddeath_date),
      death_date, 
      as.Date(NA_character_))
  ) %>%
  # discard all data before the earliest time zero for each individual
  mutate(across(c(postest_date, covidadmitted_date, coviddeath_date, death_date, noncoviddeath_date),
                ~ if_else(!is.na(.x) & start_fu_date < .x,
                          .x,
                          as.Date(NA_character_)))) 

# save outcomes data
readr::write_rds(
  data_outcomes,
  here::here("output", glue("jcvi_group_{group}"), "data", "data_outcomes.rds"),
  compress = "gz")
