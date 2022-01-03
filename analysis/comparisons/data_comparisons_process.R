################################################################################

# This script:
# reads:
## second vaccination dates, 
## data from eligible individuals based on boxes c and d
## covariates data
# derives:
## follow-up times for each comparison
## covariates updated at time_zero for each comparison
# saves:
## comparison data with covariates (one-row-per-comparison-per-individual)

################################################################################

library(tidyverse)
library(lubridate)
library(glue)

################################################################################

# create output directories ----
fs::dir_create(here::here("output", "comparisons", "data"))

# import functions ----
source(here::here("analysis", "lib", "data_process_functions.R"))

# import variable names
model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

# import data ----
# second vaccination period dates
second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds"))

# individuals eligible based on box c criteria
data_eligible_c <- readr::read_rds(
  here::here("output", "data", "data_eligible_c.rds"))

# individuals eligible based on box d criteria
data_eligible_d <- readr::read_rds(
  here::here("output", "data", "data_eligible_d.rds"))

# covariate data
data_covs <- readr::read_rds(
  here::here("output", "data", "data_covs.rds"))

################################################################################

# derive comparison arms for k comparisons ----
# define start_fu_date & end_fu_date for each comparison
comparison_arms <- function(
  b, # brand: "BNT162b2", "ChAdOx"
  a, # arm: "vax", "unvax"
  k # comparison number, k=1...n_comparisons
) {
  
  # comparison starts on d days since second vax date
  d <- 14 + (k-1)*28
  
  # function to be applied in dplyr::filter
  no_evidence_of <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
  # exclude if evidence of these variables before index_date
  exclude_if_evidence_of <- c(
    "positive_test_0_date",
    "primary_care_covid_case_0_date",
    "primary_care_suspected_covid_0_date",
    "covidadmitted_0_date",
    "coviddeath_date",
    "endoflife_date",
    "midazolam_date",
    "death_date"
  )
  
  # which split to keep for comparison k
  split_string <- if_else((k %% 2) == 0, "even", "odd") 
  
  if (a == "vax") {
    
    data <- data_eligible_c %>% 
      # keep the given brand
      filter(brand %in% b) %>%
      # start date for vax arm depends on individual's second vax date
      mutate(
        start_fu_date = covid_vax_2_date + days(d),
        end_fu_date = start_fu_date + days(28)
        ) %>%
      # remove individuals for whom end_fu_date is after study_parameters$end_date
      filter(
        end_fu_date <= study_parameters$end_date
      ) %>%
      # no third dose before start_fu_date
      filter(no_evidence_of(covid_vax_3_date, start_fu_date)) %>%
      mutate(arm = "vax") %>%
      droplevels()
      
  } else if (a == "unvax") {
    
    data <- data_eligible_d %>% 
      # keep the given brand
      filter(brand %in% b) %>%
      # start_fu_date for unvax arm depends on elig_date, region and brand
      mutate(
        start_fu_date = start_of_period + days(d),
        end_fu_date = start_fu_date + days(56)
        ) %>%
      # remove individuals for whom end_fu_date is after study_parameters$end_date
      filter(
        end_fu_date <= study_parameters$end_date
      ) %>%
      filter(
        # no first dose before start_fu_date
        no_evidence_of(covid_vax_1_date, start_fu_date),
        # only keep 50% of unvax individuals, depending on if k odd or even
        split %in% split_string
      ) %>%
      mutate(arm = "unvax") %>%
      select(-split) %>%
      droplevels()
      
  } else {
    
    stop("a must be \"vax\" or \"unvax\"")
    
  }
  
  # bind datasets from both arms and apply exclusions
  out <- data %>%
    left_join(data_covs %>%
                select(patient_id, 
                       all_of(exclude_if_evidence_of)),
              by = "patient_id") %>%
    # exclude if evidence of xxx before start_fu_date
    filter_at(
      all_of(exclude_if_evidence_of),
      all_vars(no_evidence_of(., start_fu_date))) %>%
    select(
      patient_id, jcvi_group, elig_date, region_0, brand, arm, 
      start_fu_date, end_fu_date
      ) %>% 
    mutate(comparison = k)
  
  return(out)
  
}

# check that no overlap between individuals in:
## brand arms
## vax and unvax arms
## unvax arm in odd and even comparisons
check_overlap <- function(.data) {
  
  data <- .data %>%
    select(patient_id, brand, arm, comparison) %>%
    mutate(
      name = if_else(
        as.numeric(comparison) %% 2 == 0,
        "even", 
        "odd"),
      value = TRUE) %>%
    mutate(across(name, 
                  ~ if_else(
                    arm %in% "vax",
                    as.character(arm),
                    str_c(arm, .x, sep = "_")))) %>%
    mutate(across(name, ~str_c(brand, .x, sep = "_"))) %>%
    distinct(patient_id, name, value) %>%
    pivot_wider(
      names_from = name, values_from = value
    ) %>%
    mutate(across(-patient_id, ~!is.na(.x))) %>%
    filter(rowSums(.[-1]) > 1)
  
  
  stopifnot("Overlap ids between arms or odd and even comparisons" = nrow(data) == 0)
  
  return(TRUE)
  
}

# BNT162b2 vs unvax
data_comparisons_pfizer <- bind_rows(
  lapply(
    c("vax", "unvax"),
    function(x)
      lapply(
        1:study_parameters$max_comparisons,
        function(y)
          comparison_arms(b = "BNT162b2", a = x, k = y)
      )
  )
) 

if (data_comparisons_pfizer %>% check_overlap()) {
  readr::write_rds(
    data_comparisons_pfizer,
    here::here("output", "data", "data_comparisons_pfizer.rds"),
    compress = "gz"
  )
}

# ChAdOx vs unvax
data_comparisons_az <- bind_rows(
  lapply(
    c("vax", "unvax"),
    function(x)
      lapply(
        1:study_parameters$max_comparisons,
        function(y)
          comparison_arms(b = "ChAdOx", a = x, k = y)
      )
  )
) 

if (data_comparisons_az %>% check_overlap()) {
  readr::write_rds(
    data_comparisons_az,
    here::here("output", "data", "data_comparisons_az.rds"),
    compress = "gz"
  )
}

# BNT162b2 vs ChAdOx
data_comparisons_both <- bind_rows(
  lapply(
    c("BNT162b2", "ChAdOx"),
    function(x)
      lapply(
        1:study_parameters$max_comparisons,
        function(y)
          comparison_arms(b = x, a = "vax", k = y)
      )
  )
) 

if (data_comparisons_both %>% check_overlap()) {
  readr::write_rds(
    data_comparisons_both,
    here::here("output", "data", "data_comparisons_both.rds"),
    compress = "gz"
  )
}