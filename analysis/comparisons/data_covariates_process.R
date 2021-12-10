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

# import functions ----
source(here::here("analysis", "lib", "data_process_functions.R"))

# import data ----
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_a <- readr::read_rds(
  here::here("output", "vax", "data", "data_eligible_a.rds")) %>%
  filter(jcvi_group %in% group)

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_b <- readr::read_rds(
  here::here("output", "vax", "data", "data_eligible_b.rds")) %>%
  filter(jcvi_group %in% group)

data_vax_wide <- readr::read_rds(
  here::here("output", "vax", "data", "data_wide_vax_dates.rds"))

second_vax_period_dates <- readr::read_csv(
  here::here("output", "lib", "second_vax_period_dates.csv"))

input_covs <- arrow::read_feather(
  here::here("output", "input_covs.feather")) %>%
  mutate(across(where(is.POSIXct), as.Date))

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
  select(patient_id, jcvi_group, elig_date, region_0, ethnicity, covid_vax_2_date, covid_vax_3_date, brand, start_of_period, end_of_period)

# apply eligibility criteria in box d ----
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

# derive comparison arms for k comparisons ----
# define time_zero & end_fu_date for each comparion
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
      mutate(time_zero = covid_vax_2_date + days(d)) %>%
      # no third dose before time_zero
      filter(no_evidence_of(covid_vax_3_date, time_zero)) %>%
      mutate(arm = "vax")
    
    data_unvax <- data_eligible_d %>%
      # time_zero for unvax arm depends on elig_date, region and brand
      mutate(time_zero = start_of_period + days(d)) %>%
      # no first dose before time_zero
      filter(no_evidence_of(covid_vax_1_date, time_zero)) %>%
      mutate(arm = "unvax")
    
  
  make_exclusions <- function(.data) {
    .data %>%
      left_join(input_covs %>%
                  select(patient_id, 
                         all_of(exclude_if_evidence_of), 
                         region = glue("region_{k}")),
                by = "patient_id") %>%
      # exclude if evidence of xxx before time_zero
      filter_at(all_of(exclude_if_evidence_of),
                all_vars(no_evidence_of(., time_zero))) %>%
      select(patient_id, elig_date, region, ethnicity, brand, arm, time_zero) 
  }
  
  out <- bind_rows(make_exclusions(data_vax), make_exclusions(data_unvax)) 
  
  end_fu_dates <- out %>%
    group_by(elig_date, brand, region) %>%
    summarise(end_fu_date = max(time_zero) + days(28), .groups = "keep") %>%
    ungroup() %>%
    mutate(arm = "unvax")
 
  out <- out %>%
    left_join(end_fu_dates, by = c("elig_date", "region", "brand", "arm")) %>%
    mutate(across(end_fu_date,
                  ~ if_else(arm %in% "vax", 
                            # each individual in vax arm followed up for 28 days
                            time_zero + days(28), 
                            .x)))
  
  return(out)
    
}

data_comparison_arms <- bind_rows(lapply(
  1:study_parameters$n_comparisons,
  function(x)
    comparison_arms(k=1) %>% mutate(k = x)
)) %>%
  mutate(across(arm, factor, levels = c("unvax", "vax")))


# add clinical and demographic covariates ----
strata_vars <- c("region", "elig_date", "k")
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

# derive long datasets from recurring variables ----
# imd IS NOT brand specific
imd_data <- data_comparison_arms %>%
  distinct(patient_id) %>%
  left_join(input_covs %>%
              select(patient_id, starts_with("imd")),
            by = "patient_id") %>%
  pivot_longer(cols = -patient_id) %>%
  mutate(k = as.integer(str_extract(name, "\\d+"))) %>%
  select(patient_id, k, imd = value)

# shielded IS brand specific
shielded_data <- data_comparison_arms %>%
  distinct(patient_id, k, brand, time_zero) %>%
  left_join(input_covs %>%
              select(patient_id, contains("shield")) %>%
              pivot_longer(cols = contains("shield")) %>%
              filter(!is.na(value)) %>%
              mutate(across(name, ~str_remove(.x, "_\\d+_date"))),
            by = "patient_id") %>%
  filter(!is.na(value) & value <= time_zero) %>%
  group_by(patient_id, brand, k, name) %>%
  summarise(date = max(value), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(names_from = name, values_from = date) %>%
  mutate(shielded = case_when(
    !is.na(shielded) & 
      ((!is.na(nonshielded) & (shielded > nonshielded)) |
         is.na(nonshielded)) ~ TRUE,
    TRUE ~ FALSE)) %>%
  select(patient_id, brand, k, shielded)

# bmi IS brand specific
bmi_data <- data_comparison_arms %>%
  distinct(patient_id, brand, k, time_zero) %>%
  left_join(input_covs %>%
              select(patient_id, contains("bmi")) %>%
              rename_if(is.Date,
                        ~ str_c("date_", str_extract(.x, "\\d+"))) %>%
              pivot_longer(cols =  contains(c("bmi", "date")),
                           names_sep = "_",
                           names_to = c(".value", "i")) %>%
              filter(!is.na(bmi)),
            by = "patient_id") %>%
  filter(!is.na(bmi) & date <= time_zero ) %>%
  group_by(patient_id, brand, k) %>%
  arrange(desc(date), .by_group = TRUE) %>%
  # keeps just the most recent before time_zero
  distinct(patient_id, k, .keep_all = TRUE) %>%
  ungroup() %>%
  select(patient_id, brand, k, bmi)
  
  
# add covariates to data_comparison_arms ----
data_covariates <- data_comparison_arms %>%
  # join and process covariates
  left_join(
    input_covs %>%
      select(patient_id, 
             dob, 
             all_of(demographic_vars[demographic_vars %in% names(.)]),
             all_of(ever_vars),
             all_of(clinical_vars[clinical_vars %in% names(.)]),
             all_of(end_vars)
      ), 
    by = "patient_id") %>%
  left_join(
    imd_data, 
    by = c("patient_id", "k")) %>%
  left_join(
    shielded_data, 
    by = c("patient_id", "brand", "k")) %>%
  mutate(across(
    shielded, 
    ~ if_else(is.na(.x), FALSE, .x))) %>%
  left_join(
    bmi_data, 
    by = c("patient_id", "brand", "k")) %>%
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
       .x <= time_zero ~ TRUE,
       TRUE ~ FALSE))) %>%
  rename_with(
    .fn = ~ str_remove(.x, "_date"),
    .cols = all_of(ever_vars)
    ) %>%
  mutate(
    
    age = as.numeric(time_zero - dob),
    
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
    
    efi = fct_case_when(
      is.na(efi) | (efi <= 0.12) ~ "None",
      efi <= 0.24 ~ "Mild",
      efi <= 0.36 ~ "Moderate",
      TRUE ~ "Severe"
    )
    
  ) 

readr::write_rds(
  data_covariates,
  here::here("output", glue("jcvi_group_{group}"), "data", "data_covariates.rds"),
  compress = "gz")


