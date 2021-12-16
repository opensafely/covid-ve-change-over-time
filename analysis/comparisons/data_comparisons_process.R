######################################

# This script:


######################################

# import libraries ----
library(tidyverse)
library(lubridate)
library(glue)

# import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  group <- "02"
  
} else {
  group <- args[[1]]
}

## create output directories ----
fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "data"))
fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "models"))

# import functions ----
source(here::here("analysis", "lib", "data_process_functions.R"))

model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

# import data ----
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

data_eligible_c <- readr::read_rds(
  here::here("output", "data", "data_eligible_c.rds")) %>%
  filter(jcvi_group %in% group)

data_eligible_d <- readr::read_rds(
  here::here("output", "data", "data_eligible_d.rds")) %>%
  filter(jcvi_group %in% group)

input_covs <- arrow::read_feather(
  here::here("output", "input_covs.feather")) %>%
  mutate(across(where(is.POSIXct), as.Date)) %>%
  filter(jcvi_group %in% group)

################################################################################

# derive comparison arms for k comparisons ----
# define time_zero_date & end_fu_date for each comparison
comparison_arms <- function(
  k # comparison number, k=1...K
) {
  
  # comparison starts on d days since second vax date
  d <- 14 + (k-1)*28
  
  no_evidence_of <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
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
    
    data_vax <- data_eligible_c %>%
      # start date for vax arm depends on second vax date
      mutate(time_zero_date = covid_vax_2_date + days(d)) %>%
      # no third dose before time_zero_date
      filter(no_evidence_of(covid_vax_3_date, time_zero_date)) %>%
      mutate(arm = "vax") %>%
      droplevels()
    
    data_unvax <- data_eligible_d %>%
      # time_zero_date for unvax arm depends on elig_date, region and brand
      mutate(time_zero_date = start_of_period + days(d)) %>%
      # no first dose before time_zero_date
      filter(no_evidence_of(covid_vax_1_date, time_zero_date)) %>%
      mutate(arm = "unvax") %>%
      droplevels()
    
  
  make_exclusions <- function(.data) {
    .data %>%
      left_join(input_covs %>%
                  select(patient_id, 
                         all_of(exclude_if_evidence_of), 
                         region = glue("region_{k}")),
                by = "patient_id") %>%
      # exclude if evidence of xxx before time_zero_date
      filter_at(all_of(exclude_if_evidence_of),
                all_vars(no_evidence_of(., time_zero_date))) %>%
      select(patient_id, elig_date, region, ethnicity, brand, arm, time_zero_date) 
  }
  
  out <- bind_rows(make_exclusions(data_vax), make_exclusions(data_unvax)) 
  
  # elig_date:brand:region-specific end_fu_date for unvax arm
  end_fu_dates <- out %>%
    group_by(elig_date, brand, region) %>%
    summarise(end_fu_date = max(time_zero_date) + days(28), .groups = "keep") %>%
    ungroup() %>%
    mutate(arm = "unvax")
 
  # join and derive individual-specific end_fu_date for vax arm
  out <- out %>%
    left_join(end_fu_dates, by = c("elig_date", "region", "brand", "arm")) %>%
    mutate(across(end_fu_date,
                  ~ if_else(arm %in% "vax", 
                            # each individual in vax arm followed up for 28 days
                            time_zero_date + days(28), 
                            .x))) %>% 
    mutate(comparison = k)
  
  return(out)
    
}

data_comparison_arms <- bind_rows(lapply(
  1:study_parameters$n_comparisons,
  function(x)
    comparison_arms(k=x))) %>%
  mutate(across(arm, factor, levels = c("vax", "unvax")))  %>%
  mutate(across(comparison, factor))

################################################################################
# read long datasets from recurring variables ----

# imd (index is non-brand-specific start_k)
data_imd <- input_covs %>%
  select(patient_id, starts_with("imd")) %>%
  pivot_longer(cols = -patient_id) %>%
  mutate(comparison = factor(as.integer(str_extract(name, "\\d+")))) %>%
  select(patient_id, comparison, imd = value)

# shielded (index is time_zero)
data_shielded <- data_comparison_arms %>%
  distinct(patient_id, comparison, brand, time_zero_date) %>%
  inner_join(
    bind_rows(
      readr::read_rds(
        here::here("output", "data", "data_long_shielded_dates.rds")) %>%
        select(patient_id, date) %>%
        mutate(name = "shielded"),
      readr::read_rds(
        here::here("output", "data", "data_long_nonshielded_dates.rds")) %>%
        select(patient_id, date) %>%
        mutate(name = "nonshielded")),
            by = "patient_id") %>%
  filter(date <= time_zero_date) %>%
  group_by(patient_id, brand, comparison, name) %>%
  summarise(date = max(date), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(names_from = name, values_from = date) %>%
  mutate(shielded = case_when(
    !is.na(shielded) & 
      ((!is.na(nonshielded) & (shielded > nonshielded)) |
         is.na(nonshielded)) ~ TRUE,
    TRUE ~ FALSE)) %>%
  select(patient_id, brand, comparison, shielded)

# bmi (index is time_zero)
data_bmi <- data_comparison_arms %>%
  distinct(patient_id, brand, comparison, time_zero_date) %>%
  inner_join(readr::read_rds(
    here::here("output", "data", "data_long_bmi_dates.rds")) %>%
      select(patient_id, date, bmi),
            by = "patient_id") %>%
  filter(date <= time_zero_date ) %>%
  group_by(patient_id, brand, comparison) %>%
  arrange(desc(date), .by_group = TRUE) %>%
  # keeps just the most recent before time_zero_date
  distinct(patient_id, comparison, .keep_all = TRUE) %>%
  ungroup() %>%
  select(patient_id, brand, comparison, bmi)

# hospital admissions
# TODO



################################################################################

# add clinical and demographic covariates to data_comparison_arms ----

strata_vars <- c("region", "elig_date", "comparison")
demographic_vars <- c("age", "sex", "ethnicity", "imd")
ever_vars <- c(
  "longres_date",
  "asplenia_date", 
  "bone_marrow_transplant_date", 
  "cancer_excl_lung_and_haem_date",  
  "chemo_or_radio_date",
  "chronic_cardiac_disease_date",
  "chronic_liver_disease_date",
  "current_copd_date",
  "cystic_fibrosis_date",
  "dementia_date",
  "diabetes_date",
  "dialysis_date",
  "dmards_date",
  "haematological_cancer_date",
  "heart_failure_date",
  "ld_inc_ds_and_cp_date",
  "lung_cancer_date",
  "other_heart_disease_date",
  "other_neuro_conditions_date",
  "other_respiratory_date",
  "permanant_immunosuppression_date",
  "psychosis_schiz_bipolar_date",
  "sickle_cell_disease_date",
  "solid_organ_transplantation_date",
  "temporary_immunosuppression_date"
)
clinical_vars <- c(
  "flu_vaccine",
  "efi",
  "bmi"
)
end_vars <- c(
  "coviddeath_date",
  "death_date",
  "dereg_date"
)


data_comparisons <- data_comparison_arms %>%
  # join and process covariates
  left_join(
    input_covs %>%
      select(patient_id, 
             dob, 
             all_of(demographic_vars[demographic_vars %in% names(.)]),
             all_of(ever_vars),
             all_of(clinical_vars[clinical_vars %in% names(.)]),
             all_of(end_vars)), 
    by = "patient_id") %>%
  left_join(
    data_imd, 
    by = c("patient_id", "comparison")) %>%
  mutate(across(imd, factor)) %>%
  left_join(
    data_shielded, 
    by = c("patient_id", "brand", "comparison")) %>%
  mutate(across(
    shielded, 
    ~ if_else(is.na(.x), FALSE, .x))) %>%
  left_join(
    data_bmi, 
    by = c("patient_id", "brand", "comparison")) %>%
  mutate(across(
    bmi,
    ~ fct_case_when(is.na(.x) | .x < 30 | .x >= 100 ~ "Not Obese",
                    .x < 35 ~ "Obese I (30-34.9)",
                    .x < 40 ~ "Obese II (35-39.9)",
                    TRUE ~ "Obese III (40+)"))) %>%
  mutate(across(
    all_of(ever_vars), 
     ~ case_when(
       is.na(.x) ~ FALSE,
       .x <= time_zero_date ~ TRUE,
       TRUE ~ FALSE))) %>%
  rename_with(
    .fn = ~ str_remove(.x, "_date"),
    .cols = all_of(ever_vars)
    ) %>%
  mutate(
    
    age = as.numeric(time_zero_date - dob)/365.25,
    
    any_immunosuppression = (permanant_immunosuppression | 
                               asplenia | 
                               dmards | 
                               solid_organ_transplantation |
                               sickle_cell_disease | 
                               temporary_immunosuppression | 
                               bone_marrow_transplant | 
                               chemo_or_radio),
    
    multimorb =
      (bmi %in% c("Obese II (35-39.9)", "Obese III (40+)")) +
      (chronic_cardiac_disease | heart_failure | other_heart_disease) +
      (dialysis) +
      (diabetes) +
      (chronic_liver_disease)+
      (current_copd | other_respiratory)+
      (lung_cancer | haematological_cancer | cancer_excl_lung_and_haem)+
      (any_immunosuppression)+
      (dementia | other_neuro_conditions)+
      (ld_inc_ds_and_cp)+
      (psychosis_schiz_bipolar),
    multimorb = cut(
      multimorb, 
      breaks = c(0, 1, 2, 3, 4, Inf),
      labels=c("0", "1", "2", "3", "4+"), 
      right=FALSE),
    
    flu_vaccine = flu_vaccine == 1,
    
    efi = fct_case_when(
      is.na(efi) | (efi <= 0.12) ~ "None",
      efi <= 0.24 ~ "Mild",
      efi <= 0.36 ~ "Moderate",
      TRUE ~ "Severe"
    )
    
  ) %>%
  select(
   patient_id, elig_date, region, ethnicity, brand, arm, 
   time_zero_date, end_fu_date, comparison,
   dob, 
   unname(unlist(model_varlist))
  ) %>%
  droplevels()

readr::write_rds(
  data_comparisons,
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds"),
  compress = "gz")
