################################################################################

# This script:
# - reads the processed data
# - applies eligibility criteria from boxes a and b of Figure 3 in protocol
# - saves preocessed data from eligible individuals

################################################################################

# import libraries ----
library(tidyverse)
library(lubridate)
library(glue)

study_parameters <- readr::read_rds(here::here("output", "lib", "study_parameters.rds"))

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_a <- readr::read_rds(
  here::here("output", "data", "data_eligible_a.rds")) 

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_b <- readr::read_rds(
  here::here("output","data", "data_eligible_b.rds")) 

# read wide vaccine dates data
data_vax_wide <- readr::read_rds(
  here::here("output", "data", "data_wide_vax_dates.rds"))

# read second vax period dates and filter to brands with enough individuals
second_vax_period_dates <- readr::read_csv(
  here::here("output", "lib", "second_vax_period_dates.csv")) %>%
  filter(n_in_period > study_parameters$n_threshold)

################################################################################
# apply eligibility criteria in box c ----
data_eligible_c <- data_eligible_b %>%
  left_join(data_vax_wide, 
            by = "patient_id") %>%
  mutate(brand = case_when(covid_vax_2_brand %in% "pfizer" ~ "BNT162b2",
                           covid_vax_2_brand %in% "az" ~ "BNT162b2",
                           TRUE ~ NA_character_)) %>%
  select(-ends_with("_brand")) %>%
  left_join(second_vax_period_dates, 
            by = c("jcvi_group", "elig_date", "region_0", "brand")) %>%
  filter(
    # second dose during second vax period
    start_of_period <= covid_vax_2_date,
    covid_vax_2_date <= end_of_period,
    # enough individuals in second vax period for a given elig_date and brand to include comparison (i.e. summed over regions)
    n_in_period >= 100) %>%
  select(patient_id, jcvi_group, elig_date, region_0, ethnicity, 
         covid_vax_2_date, covid_vax_3_date, brand, 
         start_of_period, end_of_period) %>%
  droplevels()

readr::write_rds(
  data_eligible_c,
  here::here("output", "data", "data_eligible_c.rds"),
  compress = "gz")

################################################################################
# apply eligibility criteria in box d ----
data_eligible_d <- data_eligible_a %>%
  left_join(data_vax_wide %>%
              select(-ends_with("_brand")),
            by = "patient_id") %>%
  # creates 1 row per brand, so some individuals in the unvax arm will have 2 rows 
  left_join(second_vax_period_dates, 
            by = c("jcvi_group", "elig_date", "region_0")) %>%
  # remove individuals who had received any vaccination before the start of the second vax period
  filter(
    is.na(covid_vax_1_date) | covid_vax_1_date >= start_of_period
  ) %>%
  droplevels()

readr::write_rds(
  data_eligible_d,
  here::here("output", "data", "data_eligible_d.rds"),
  compress = "gz")
