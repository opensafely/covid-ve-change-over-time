################################################################################

# This script:
# - processes the outcome data

################################################################################

library(tidyverse)
library(lubridate)
library(glue)


################################################################################
# read datasets ----
data_long_postest_dates <- readr::read_rds(
  here::here("output", "data", "data_long_postest_dates.rds"))  %>%
  select(patient_id, postest_date = date)

data_long_covidadmitted_dates <- readr::read_rds(
  here::here("output", "data", "data_long_covidadmitted_dates.rds")) %>%
  select(patient_id, covidadmitted_date = date)

data_covs <- readr::read_rds(
  here::here("output", "data", "data_covs.rds")) %>%
  select(patient_id, coviddeath_date, death_date, dereg_date) %>%
  mutate(across(ends_with("date"), as.Date))

################################################################################
# process outcomes data ----
for (arm in c("unvax", "BNT162b2", "ChAdOx")) {
  
  cat(glue("----process outcomes data for arm: {arm}----"), "\n")
  
  # read comparisons data
  # keep only patient_id and earliest start_fu_date 
  # (i.e. from comparison 1 for vax and odd unvax, comparison 2 for even unvax)
  data_ids <- readr::read_rds(
    here::here("output", "data", glue("data_comparisons_{arm}.rds"))) %>%
    arrange(patient_id, start_fu_date) %>%
    select(patient_id, start_fu_date) %>%
    distinct(patient_id, .keep_all = TRUE)
  
  ## join outcomes data
  data_outcomes <- data_ids %>%
    left_join(
      data_long_postest_dates,
      by = "patient_id"
    ) %>%
    left_join(
      data_long_covidadmitted_dates,
      by = "patient_id"
    ) %>%
    left_join(
      data_covs,
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
    # discard all data before start_fu_date for each individual
    mutate(across(c(postest_date, covidadmitted_date, coviddeath_date, death_date, noncoviddeath_date),
                  ~ if_else(!is.na(.x) & start_fu_date < .x,
                            .x,
                            as.Date(NA_character_)))) 
  
  # save outcomes data
  readr::write_rds(
    data_outcomes,
    here::here("output", "data", glue("data_outcomes_{arm}.rds")),
    compress = "gz")
  
}
