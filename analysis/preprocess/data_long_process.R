###############################################################################

# What this script does:

# reads covariate data from eligible individuals
# create long (one-row-per-event) datasets for recurring variables
# saves long datasets

###############################################################################
library(tidyverse)
library(glue)

## preliminaries
# import covariates data
data_covs <- readr::read_rds(
  here::here("output", "data", "data_covs.rds"))

# import functions ----
source(here::here("analysis", "lib", "data_process_functions.R"))

###############################################################################
## create one-row-per-event datasets

# shielded
data_pr_shielded <- data_covs %>%
  select(patient_id,
         matches("^shielded\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "shielded_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

readr::write_rds(
  data_pr_shielded, 
  here::here("output", "data", "data_long_shielded_dates.rds"), 
  compress="gz")

###############################################################################
# nonshielded
data_pr_nonshielded <- data_covs %>%
  select(patient_id,
         matches("^nonshielded\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "nonshielded_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

readr::write_rds(
  data_pr_nonshielded, 
  here::here("output", "data", "data_long_nonshielded_dates.rds"), 
  compress="gz")

###############################################################################
# bmi
data_pr_bmi <- data_covs %>%
  select(patient_id,
         matches("^bmi\\_\\d+")) %>%
  rename_at(vars(contains("date")),
            ~ str_c("date_", str_extract(.x, "\\d+"))) %>%
  pivot_longer(
    cols = -patient_id,
    names_sep = "_",
    names_to = c(".value", "bmi_index"),
    values_drop_na = TRUE) %>%
  mutate(bmi = fct_case_when(
    bmi < 30 | bmi >=100 ~ "Not obese", # this cat includes missing and clinically implausible values
    bmi >= 30 & bmi < 35 ~ "Obese I (30-34.9)",
    bmi >= 35 & bmi < 40 ~ "Obese II (35-39.9)",
    bmi >= 40 & bmi < 100 ~ "Obese III (40+)",
    TRUE ~ NA_character_
  )) %>%
  arrange(patient_id, date) 

readr::write_rds(
  data_pr_bmi, 
  here::here("output", "data", "data_long_bmi_dates.rds"), 
  compress="gz")

###############################################################################
# suspected covid
data_pr_suspected_covid <- data_covs %>%
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

readr::write_rds(
  data_pr_suspected_covid, 
  here::here("output", "data", "data_long_pr_suspected_covid_dates.rds"), 
  compress="gz")

###############################################################################
# probable covid
data_pr_probable_covid <- data_covs %>%
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

readr::write_rds(
  data_pr_probable_covid, 
  here::here("output", "data", "data_long_pr_probable_covid_dates.rds"),
  compress="gz")

###############################################################################
# positive test
data_postest <- data_covs %>%
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
  data_postest, 
  here::here("output", "data", "data_long_postest_dates.rds"), 
  compress="gz")

###############################################################################
# covid admission
data_covidadmitted <- data_covs %>%
  select(patient_id, 
         matches("^covidadmitted\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "covidadmitted_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

readr::write_rds(
  data_covidadmitted, 
  here::here("output", "data", "data_long_covidadmitted_dates.rds"), 
  compress="gz")

# ###############################################################################
# # hospital admissions
# data_admissions <- data_covs %>%
#   select(patient_id, 
#          matches("^admitted\\_unplanned\\_\\d+\\_date"), 
#          matches("^discharged\\_unplanned\\_\\d+\\_date")) %>%
#   pivot_longer(
#     cols = -patient_id,
#     names_to = c(".value", "index"),
#     names_pattern = "^(.*)_(\\d+)_date",
#     values_drop_na = TRUE
#   ) %>%
#   select(patient_id, index, 
#          admitted_date=admitted_unplanned, 
#          discharged_date = discharged_unplanned) %>%
#   arrange(patient_id, admitted_date)
# 
# readr::write_rds(
#   data_admissions, 
#   here::here("output", "data", "data_long_admission_dates.rds"), 
#   compress="gz")
# 
# ###############################################################################
# # infectious hospital admissions
# data_admissions_infectious <- data_covs %>%
#   select(patient_id, 
#          matches("^admitted\\_unplanned\\_infectious\\_\\d+\\_date"),
#          matches("^discharged\\_unplanned\\_infectious\\_\\d+\\_date")) %>%
#   pivot_longer(
#     cols = -patient_id,
#     names_to = c(".value", "index"),
#     names_pattern = "^(.*)_(\\d+)_date",
#     values_drop_na = TRUE
#   ) %>%
#   select(patient_id, index, 
#          admitted_date=admitted_unplanned_infectious, 
#          discharged_date = discharged_unplanned_infectious) %>%
#   arrange(patient_id, admitted_date)
# 
# readr::write_rds(
#   data_admissions_infectious, 
#   here::here("output", "data", "data_long_admission_infectious_dates.rds"), 
#   compress="gz")
# 
# ###############################################################################
# # noninfectious hospital admissions
# #remove infectious admissions from all admissions data
# data_admissions_noninfectious <- anti_join(
#   data_admissions,
#   data_admissions_infectious,
#   by = c("patient_id", "admitted_date", "discharged_date")
# )
# 
# readr::write_rds(
#   data_admissions_noninfectious, 
#   here::here("output", "data", "data_long_admission_noninfectious_dates.rds"), 
#   compress="gz")




