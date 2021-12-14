library(tidyverse)
library(lubridate)
library(glue)

source(here::here("analysis", "lib", "dummy_data_functions.R"))

set.seed(5476)

# the variables generated in grouping_variables.py are same as in dummy_data_vax
dummy_data_vax <- arrow::read_feather(
  here::here("analysis", "vax", "dummy_data_vax.feather"))

study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

# start and end dates
start_dates <- readr::read_csv(here::here("output", "lib", "start_dates.csv"))
end_dates <- readr::read_csv(here::here("output", "lib", "end_dates.csv"))

translate_to_R <- function(.data) {
  .data %>%
    mutate(across(condition, ~ str_replace_all(.x, "elig_date = ", "elig_date ==\'"))) %>%
    mutate(across(condition, ~ str_replace_all(.x, "region_0 = ", "region_0 =="))) %>%
    mutate(across(condition, ~ str_replace_all(.x, " AND", "\' &"))) %>%
    mutate(across(condition, ~ str_replace_all(.x, "OR", "|"))) %>%
    mutate(across(condition, ~ str_replace(.x, "DEFAULT", "TRUE"))) 
}
# translate conditions to R
start_dates <- start_dates %>% translate_to_R()
end_dates <- end_dates %>% translate_to_R()

# date vars (incidences are higher than reality to ensure code runs on smaller sample)
date_vars_ever <- c("chronic_cardiac_disease_date",
                    "heart_failure_date",
                    "other_heart_disease_date",
                    "diabetes_date",
                    "dialysis_date", 
                    "chronic_liver_disease_date", 
                    "current_copd_date",
                    "ld_inc_ds_and_cp_date", 
                    "cystic_fibrosis_date", 
                    "other_respiratory_date",
                    "lung_cancer_date",
                    "haematological_cancer_date",
                    "cancer_excl_lung_and_haem_date", 
                    "chemo_or_radio_date",
                    "solid_organ_transplantation_date",
                    "bone_marrow_transplant_date",
                    "sickle_cell_disease_date", 
                    "permanant_immunosuppression_date",
                    "temporary_immunosuppression_date", 
                    "asplenia_date", 
                    "dmards_date", 
                    "dementia_date", 
                    "other_neuro_conditions_date", 
                    "psychosis_schiz_bipolar_date")

date_vars_recent <- c("positive_test_0_date", 
                      "primary_care_covid_case_0_date", 
                      "primary_care_suspected_covid_0_date",
                      "covidadmitted_0_date",
                      "death_date",
                      "longres_date", 
                      "endoflife_date", 
                      "midazolam_date",
                      "coviddeath_date", 
                      "dereg_date")


dummy_data_covs <- dummy_data_vax %>%
  select(patient_id, age_1, age_2, sex, jcvi_group, elig_date, region_0) %>%
  mutate(across(ends_with("_date"), as.Date)) %>%
  # start and end dates for second vax period
  var_category(
    start_1_date,
    categories = start_dates$start_1_date,
    conditions = start_dates$condition
  ) %>%
  var_category(
    end_1_date,
    categories = end_dates$end_1_date,
    conditions = end_dates$condition
  ) %>%
  var_binary(name = flu_vaccine, incidence = 0.3) %>%
  # date vars  ever
  bind_cols(
    pmap(
      list(a = date_vars_ever, 
           b = rep(0.2, length(date_vars_ever))),
      function(a,b) 
        var_date(
          .data = ., 
          name = !! a,
          incidence = b,
          keep_vars = FALSE
        ))) %>%
  bind_cols(
    pmap(
      list(a = date_vars_recent, 
           b = rep(0.2, length(date_vars_recent))),
           function(a,b) 
             var_date(
               .data = ., 
               name = !! a,
               incidence = b,
               earliest="2020-11-01",
               latest="2021-12-31",
               keep_vars = FALSE
             ))) %>%
  mutate(across(death_date, 
                ~if_else(
                  !is.na(coviddeath_date), 
                  coviddeath_date,
                  NA_Date_))) %>%
  bind_cols(vars_bmi_recurrent(.data = .)) %>%
  bind_cols(vars_region_and_imd_recurrent(.data = ., K=K)) %>%
  bind_cols(
    var_date_recurrent(
      .data = ., 
      name_string = "shielded", 
      incidence = 0.2,
      r = study_parameters$recur_shielded)) %>%
  bind_cols(
    var_date_recurrent(
      .data = .,
      name_string = "nonshielded", 
      incidence = 0.1,
      r = study_parameters$recur_shielded)) %>%
  bind_cols(
    var_date_recurrent(
      .data = ., 
      name_string = "admitted_unplanned", 
      incidence = 0.1,
      r = study_parameters$recur_admissions)) %>%
  bind_cols(
    var_date_recurrent(
      .data = ., 
      name_string = "admitted_unplanned_infectious", 
      incidence = 0.1,
      r = study_parameters$recur_admissions)) %>%
  # all 1st of the month as dob YYYY-MM in OpenSAFELY
  mutate(
    dob = as.POSIXct(sample(
      x = seq(as.Date("1930-01-01"), as.Date("2006-01-01"), by = "1 month"), 
      size = nrow(.),
      replace = TRUE))) %>%
  mutate(
    efi = if_else(
      rbernoulli(n = nrow(.), p = 0.99),
      rnorm(n = nrow(.), mean = 0.2, sd = 0.09),
      NA_real_)) %>%
  mutate(across(contains("_date"), as.POSIXct)) 

dummy_data_covs <- dummy_data_covs %>%
  bind_cols(
    lapply(
      0:study_parameters$recur_admissions,
      function(x)
        dummy_data_covs %>%
        transmute(
          !! glue("discharged_unplanned_{x}_date") := 
            !! sym(glue("admitted_unplanned_{x}_date")) +
            sample(x=1:50, size = nrow(.), replace = TRUE),
          !! glue("discharged_unplanned_infectious_{x}_date") := 
            !! sym(glue("admitted_unplanned_infectious_{x}_date")) + 
            sample(x=1:50, size = nrow(.), replace = TRUE))
    )) %>%
  droplevels()


arrow::write_feather(dummy_data_covs, here::here("analysis", "covs", "dummy_data_covs.feather"))
