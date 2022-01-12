################################################################################

# This script:
# - reads the processed data
# - applies eligibility criteria from boxes c and d of Figure 3 in protocol
# - saves processed data from eligible individuals

################################################################################

# import libraries ----
library(tidyverse)
library(lubridate)
library(glue)

# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

# individuals who are eligible based on criteria in box a of Figure 3 on protocol
data_eligible_a <- readr::read_rds(
  here::here("output", "data", "data_eligible_a.rds")) 

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_b <- readr::read_rds(
  here::here("output","data", "data_eligible_b.rds")) 

# read wide vaccine dates data
data_vax_wide <- readr::read_rds(
  here::here("output", "data", "data_wide_vax_dates.rds"))

# read second vax period dates 
second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds"))

# covariate data
data_covs <- readr::read_rds(
  here::here("output", "data", "data_covs.rds")) %>%
  # date of first evidence of covid
  left_join(
    readr::read_rds(
      here::here("output", "data", "data_covid_any.rds")),
    by = "patient_id"
  ) %>%
  select(patient_id, endoflife_date, midazolam_date, covid_any_date)

################################################################################
# apply eligibility criteria in box c ----
data_eligible_c <- data_eligible_b %>%
  left_join(data_vax_wide, 
            by = "patient_id") %>%
  # keep brand of interest 
  # (already applied condition that 1st and 2nd doses are the same)
  mutate(brand = case_when(covid_vax_2_brand %in% "pfizer" ~ "BNT162b2",
                           covid_vax_2_brand %in% "az" ~ "ChAdOx",
                           TRUE ~ NA_character_)) %>%
  select(-ends_with("_brand")) %>%
  # right join to keep only the jcvi_group:elig_date:region_0:brands
  # with > n_threshold individuals vaccinated during the period
  right_join(second_vax_period_dates, 
            by = c("jcvi_group", "elig_date", "region_0", "brand")) %>%
  filter(
    # second dose during second vax period
    start_of_period <= covid_vax_2_date,
    covid_vax_2_date <= end_of_period) %>%
  select(patient_id, jcvi_group, elig_date, region_0, ethnicity, 
         covid_vax_2_date, covid_vax_3_date, brand, 
         start_of_period, end_of_period) %>%
  droplevels()

################################################################################
# apply eligibility criteria in box d ----

# set seed so that 50:50 split reproducible
set.seed(study_parameters$seed)

data_eligible_d <- data_eligible_a %>%
  # randomly split the unvax 50:50 
  # one group for odd comparisons, one for even
  # so that no overlap in follow-up time across comparisons
  mutate(split = factor(
    rbernoulli(nrow(.), p=0.5), 
    labels = c("odd", "even"))) %>%
  left_join(data_vax_wide %>%
              select(-ends_with("_brand")),
            by = "patient_id") %>%
  # right join to keep only the jcvi_group:elig_date:region_0:brands
  # with > n_threshold individuals vaccinated during the period
  # creates 1 row per brand, so some individuals in the unvax arm will have 2 rows 
  right_join(second_vax_period_dates, 
             by = c("jcvi_group", "elig_date", "region_0")) %>%
  # remove individuals who had received any vaccination before the start of the second vax period
  filter(
    is.na(covid_vax_1_date) | covid_vax_1_date >= start_of_period
  ) %>%
  select(patient_id, jcvi_group, elig_date, region_0, ethnicity, 
         covid_vax_1_date, brand, 
         start_of_period, end_of_period, split) %>%
  droplevels()

################################################################################
# apply eligibility criteria in box e ----

exclusion_e <- function(.data) {
  
  # function to be applied in dplyr::filter
  no_evidence_of <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
  .data %>%
    left_join(data_covs, by = "patient_id") %>%
    filter(
      no_evidence_of(endoflife_date, start_of_period),
      no_evidence_of(midazolam_date, start_of_period),
      no_evidence_of(covid_any_date, start_of_period - weeks(2))
    ) %>%
    select(-all_of(names(data_covs)[!names(data_covs) %in% "patient_id"]))
    
}

data_eligible_e_vax <- data_eligible_c %>% exclusion_e()
data_eligible_e_unvax <- data_eligible_d %>% exclusion_e()

readr::write_rds(
  data_eligible_e_vax,
  here::here("output", "data", "data_eligible_e_vax.rds"),
  compress = "gz")

readr::write_rds(
  data_eligible_e_unvax,
  here::here("output", "data", "data_eligible_e_unvax.rds"),
  compress = "gz")
