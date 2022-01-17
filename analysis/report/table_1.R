library(tidyverse)
library(glue)
library(lubridate)

################################################################################

fs::dir_create(here::here("output", "report", "tables"))

## import study_parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))


# individuals eligible based on box c criteria
data_eligible_e_vax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_vax.rds"))  %>%
  mutate(arm = brand) %>%
  select(patient_id, start_of_period, arm)

# individuals eligible based on box d criteria
data_eligible_e_unvax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_unvax.rds"))  %>%
  mutate(arm = "unvax") %>%
  select(patient_id, start_of_period, arm)


model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds"))

# processed data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))

# read script for processing covariates for comparison 1
source(here::here("analysis", "lib", "process_covariates.R"))

################################################################################
# function to be applied in dplyr::filter
no_evidence_of <- function(cov_date, index_date) {
  is.na(cov_date) | index_date < cov_date
}

censor_vars <- c("death_date", "dereg_date")

data_comparison_1 <- bind_rows(
  data_eligible_e_vax,
  data_eligible_e_unvax
  ) %>%
  left_join(data_processed,
            by = "patient_id") %>%
  mutate(
    # start date of comparison 1 
    start_fu_date = start_of_period + days(14),
    end_fu_date = start_fu_date + days(28)
  ) %>%
  select(-start_of_period) %>%
  # remove if death or dereg before start_of_period
  filter_at(
    all_of(censor_vars),
    all_vars(no_evidence_of(., start_fu_date))) %>%
  ungroup() %>%
  select(patient_id, jcvi_group, elig_date, region_0, arm,
         start_fu_date, end_fu_date) %>%
  mutate(comparison = 1) 

################################################################################

data_tables <- data_comparison_1 %>%
  process_covariates() %>%
  select(patient_id, arm, region, jcvi_group, subgroup,
         starts_with(unname(unlist(model_varlist)))) %>% 
  group_split(subgroup)

################################################################################

summary_var <- function(.data, var) {
  out <- .data %>%
    group_by(arm, !! sym(var)) %>%
    count() %>%
    ungroup(!! sym(var)) %>%
    mutate(arm_total = sum(n)) %>%
    ungroup() %>%
    mutate(percent = round(100*n/arm_total,0)) %>%
    mutate(across(n, ~scales::comma(.x, accuracy = 1))) %>%
    mutate(value = glue("{n} ({percent}%)")) %>%
    select(arm, !! sym(var), value) %>%
    pivot_wider(
      names_from = arm, 
      values_from = value
    ) %>%
    mutate(variable = var) %>%
    rename("category" = var)
  
  if (is.logical(out$category)) {
    out <- out %>%
      filter(category) %>%
      mutate(across(category, ~ "yes"))
  }
  
  return(out)
  
}

summary_vars <- names(data_tables[[1]] %>% select(-patient_id, -arm))

bind_rows(lapply(
  summary_vars,
  function(x)
    data_tables[[1]] %>% summary_var(var = x)
))

