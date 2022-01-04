library(tidyverse)
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
for (b in c("BNT162b2", "ChAdOx")) {
  for (a in c("vax", "unvax")) {
    
    cat(glue("----process outcomes data for {b} {a}----"), "\n")
    
    # read comparisons data
    # keep only patient_id and start_fu_date for comparison 1
    data_ids <- readr::read_rds(
      here::here("output", "data", glue("data_comparisons_{b}_{a}.rds"))) %>%
      filter(comparison %in% "1") %>%
      distinct(patient_id, start_fu_date) 
    
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
    
    
    
    # print the number of samples for which upstream outcomes missing
    # e.g. if postest missing but covidadmitted nonmissing
    data_outcomes %>%
      select(patient_id, postest_date, covidadmitted_date, coviddeath_date) %>%
      mutate(across(-patient_id, ~!is.na(.x))) %>%
      group_by(postest_date, covidadmitted_date, coviddeath_date) %>% 
      count() %>%
      ungroup() %>%
      # round n to closest 10
      mutate(n = round(n, -1)) %>%
      print(n=Inf)
    
    
    # combine outcomes where upstream outcome missing 
    # impute using median gap between outcomes
    median_times_between_outcomes <- data_outcomes %>%
      transmute(
        postest_covidadmitted = as.integer(covidadmitted_date - postest_date),
        postest_coviddeath = as.integer(coviddeath_date - postest_date),
        covidadmitted_coviddeath = as.integer(coviddeath_date - covidadmitted_date)
      ) %>%
      summarise(across(everything(), ~median(.x, na.rm = TRUE))) 
    
    
    data_outcomes_combined <- data_outcomes %>%
      mutate(across(postest_date,
                    ~ case_when(
                      !is.na(.x) ~ .x,
                      !is.na(covidadmitted_date) ~ covidadmitted_date - median_times_between_outcomes$postest_covidadmitted,
                      !is.na(coviddeath_date) ~ coviddeath_date - median_times_between_outcomes$postest_coviddeath,
                      TRUE ~ .x            
                    ))) %>%
      mutate(across(covidadmitted_date,
                    ~ case_when(
                      !is.na(.x) ~ .x,
                      !is.na(coviddeath_date) ~ coviddeath_date - median_times_between_outcomes$covidadmitted_coviddeath,
                      TRUE ~ .x            
                    )))
    
    # save outcomes data
    readr::write_rds(
      data_outcomes_combined,
      here::here("output", "data", glue("data_outcomes_{b}_{a}.rds")),
      compress = "gz")
    
    
  }
}
