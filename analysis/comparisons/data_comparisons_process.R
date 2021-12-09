######################################

# This script:


######################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_a <- readr::read_rds(
  here::here("output", "vax", "data", "data_eligible_a.rds"))

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_b <- readr::read_rds(
  here::here("output", "vax", "data", "data_eligible_b.rds"))

data_vax_wide <- readr::read_rds(
  here::here("output", "vax", "data", "data_wide_vax_dates.rds"))

second_vax_period_dates <- readr::read_csv(
  here::here("output", "lib", "second_vax_period_dates.csv"))

dummy_data_covs <- arrow::read_feather(
  here::here("analysis", "covs", "dummy_data_covs.feather"))

# apply eligiblity criteria in box c
data_eligible_c <- data_eligible_b %>%
  left_join(data_vax_wide, 
            by = "patient_id") %>%
  mutate(brand = case_when(covid_vax_2_brand %in% "pfizer" ~ "BNT162b2",
                           covid_vax_2_brand %in% "az" ~ "BNT162b2",
                           TRUE ~ NA_character_)) %>%
  select(-ends_with("_brand")) %>%
  left_join(second_vax_period_dates, 
            by = c("elig_date", "region_0", "brand")) %>%
  filter(
    # second dose during second vax period
    start_of_period <= covid_vax_2_date,
    covid_vax_2_date <= end_of_period,
    # enough individuals in second vax period for a given elig_date and brand to include comparison (i.e. summed over regions)
    n_in_period >= 100) %>%
  select(patient_id, jcvi_group, elig_date, region_0, covid_vax_2_date, covid_vax_3_date, brand, start_of_period, end_of_period)

# apply eligibility criteria in box d
data_eligible_d <- data_eligible_a %>%
  left_join(data_vax_wide %>%
              select(-ends_with("_brand")),
            by = "patient_id") %>%
  # creates 2 rows per individual, 1 for each brand
  left_join(second_vax_period_dates, 
            by = c("elig_date", "region_0")) %>%
  # remove individuals who had received any vaccination before the start of the second vax period
  filter(
    is.na(covid_vax_1_date) | covid_vax_1_date >= start_of_period
  )

######

input_covs <- arrow::read_feather(
  here::here("output", "input_covs.feather"))

comparison_arms <- function(
  j, # jcvi_group
  b, # vaccine brand
  k # comparison number, k=1...K
) {
  
  # comparison starts on d days since second vax date
  d <- 14 + (k-1)*28
  
  no_evidence_of <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
  exclude_if_evidence_of <- c(
    "positive_test_date",
    "primary_care_covid_case_date",
    "primary_care_suspected_covid_date",
    "covidadmitted_date",
    "coviddeath_date",
    "endoflife_date",
    "midazolam_date",
    "death_date"
  )
    
    data_vax <- data_eligible_c %>%
      filter(
        brand %in% b,
        jcvi_group %in% j
        ) %>%
      # start date for vax arm depends on second vax date
      mutate(time_zero = covid_vax_2_date + days(d)) %>%
      # no third dose before time_zero
      filter(no_evidence_of(covid_vax_3_date, time_zero)) %>%
      mutate(arm = "vax")
    
    data_unvax <- data_eligible_d %>%
      filter(
        brand %in% b,
        jcvi_group %in% j
      ) %>%
      # time_zero for unvax arm depends on elig_date, region and brand
      mutate(time_zero = start_of_period + days(d)) %>%
      # no first dose before time_zero
      filter(no_evidence_of(covid_vax_1_date, time_zero)) %>%
      mutate(arm = "unvax")
    
  
  make_exclusions <- function(.data) {
    .data %>%
      left_join(input_covs %>%
                  select(patient_id, all_of(exclude_if_evidence_of)),
                by = "patient_id") %>%
      # exclude if evidence of xxx before time_zero
      filter_at(all_of(exclude_if_evidence_of),
                all_vars(no_evidence_of(., time_zero))) %>%
      select(patient_id, arm, time_zero) 
  }
  
  out <- bind_rows(make_exclusions(data_vax), make_exclusions(data_unvax)) 
  
  end_fu_date <- max(out[out$arm == "vax", ]$time_zero)
  
  out <- out %>%
    mutate(end_fu_date = if_else(arm %in% "vax",
                              # each individual in vax arm followed up for 28 days
                              time_zero + days(28),
                              # each individual in unvax arm followed up until last end date in vax arm
                              end_fu_date + days(28)))
  
  return(out)
    
}

comparison_arms(j = "02", b = "BNT162b2", k=1)

