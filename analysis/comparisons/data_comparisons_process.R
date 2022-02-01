################################################################################

# This script:
# - generates 1 row per individual per comparison
# - derives comparison start and end dates
# - derives covariates

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
# individuals eligible based on box c criteria
data_eligible_e_vax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_vax.rds"))

# individuals eligible based on box d criteria
data_eligible_e_unvax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_unvax.rds"))

# processed data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))
outcomes <- outcomes[outcomes != "anytest"]

################################################################################

censor_vars <- c(
  "death_date",
  "dereg_date"
)

# derive comparison arms for k comparisons ----
# define start_fu_date & end_fu_date for each comparison
comparison_arms <- function(
  arm, # arm: "BNT162b2", "ChAdOx" "unvax"
  K = study_parameters$max_comparisons 
) {
  
  if (arm %in% "ChAdOx") {
    jcvi_groups_keep <- str_pad(as.character(2:10), width=2, side="left", pad="0")
  } else {
    jcvi_groups_keep <- str_pad(as.character(2:12), width=2, side="left", pad="0")
  }

  
  comparison_k <- function(k) {
    
    # comparison starts on d days since second vax date
    d <- 14 + (k-1)*28
    
    # function to be applied in dplyr::filter
    no_evidence_of <- function(cov_date, index_date) {
      is.na(cov_date) | index_date < cov_date
    }
    
    # which split to keep for comparison k
    split_string <- if_else((k %% 2) == 0, "even", "odd") 
    
    if (arm %in% c("BNT162b2", "ChAdOx")) {
      
      data <- data_eligible_e_vax %>% 
        # keep the given brand
        filter(brand %in% arm) %>%
        rename("arm" = "brand") %>%
        # start date for vax arm depends on individual's second vax date
        mutate(
          start_fu_date = covid_vax_2_date + days(d),
          # censor at study_parameters$end_date
          end_fu_date = pmin(start_fu_date + days(28), as.Date(study_parameters$end_date))
        ) %>%
        # remove individuals for whom start_fu_date is on or after study_parameters$end_date
        filter(
          start_fu_date < as.Date(study_parameters$end_date)
        ) %>%
        # no third dose before start_fu_date
        filter(no_evidence_of(covid_vax_3_date, start_fu_date)) %>%
        droplevels()
      
    } else if (arm %in% "unvax") {
      
      data <- data_eligible_e_unvax %>% 
        # start_fu_date for unvax arm depends on elig_date, region and brand
        mutate(
          start_fu_date = start_of_period + days(d),
          end_fu_date = pmin(start_fu_date + days(56), as.Date(study_parameters$end_date))
        ) %>%
        # remove individuals for whom start_fu_date is on or after study_parameters$end_date
        filter(
          start_fu_date < as.Date(study_parameters$end_date)
        ) %>%
        filter(
          # no first dose before start_fu_date
          no_evidence_of(covid_vax_1_date, start_fu_date),
          # only keep 50% of unvax individuals, depending on if k odd or even
          split %in% split_string
        ) %>%
        mutate(arm = arm) %>%
        select(-split) %>%
        droplevels()
      
    } else {
      
      stop("arm must be \"unvax\", \"BNT162b2\", \"ChAdOx\"")
      
    }
    
    # apply exclusions
    out_k <- data %>%
      filter(jcvi_group %in% jcvi_groups_keep) %>%
      left_join(data_processed %>%
                  select(patient_id, 
                         all_of(censor_vars)),
                by = "patient_id") %>%
      # exclude if evidence of xxx before start_fu_date
      filter_at(
        all_of(censor_vars),
        all_vars(no_evidence_of(., start_fu_date))) %>%
      select(
        patient_id, jcvi_group, elig_date, region, arm, 
        start_fu_date, end_fu_date
      ) %>% 
      mutate(comparison = k)
    
    return(out_k)
    
  }
  
  # apply comparison_k for 1:K
  out <- bind_rows(lapply(
    1:K, 
    comparison_k
  )) %>%
    mutate(across(comparison, factor, levels = 1:K))
  
  return(out)
  
}

# process covariates ----
# read script for processing covariates for each comparison
source(here::here("analysis", "lib", "process_covariates.R"))

################################################################################
# generate and save datasets

for (arm in c("unvax", "BNT162b2", "ChAdOx")) {
  
    data_comparisons <- comparison_arms(arm = arm) %>%
      process_covariates()
    
    readr::write_rds(
      data_comparisons,
      here::here("output", "comparisons", "data", glue("data_comparisons_{arm}.rds")),
      compress = "gz"
    )
    
}

