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

# import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  group <- "02"
} else {
  group <- args[[1]]
}

################################################################################

# create output directories ----
fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "data"))

# import functions ----
source(here::here("analysis", "lib", "data_process_functions.R"))

# import variable names
model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

# import data ----
# second vaccination period dates
second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds")) %>%
  filter(jcvi_group %in% group, include) 

# individuals eligible based on box c criteria
data_eligible_c <- readr::read_rds(
  here::here("output", "data", "data_eligible_c.rds")) %>%
  filter(jcvi_group %in% group)

# individuals eligible based on box d criteria
data_eligible_d <- readr::read_rds(
  here::here("output", "data", "data_eligible_d.rds")) %>%
  filter(jcvi_group %in% group)

# covariate data
data_covs <- readr::read_rds(
  here::here("output", "data", "data_covs.rds"))

################################################################################
# for the given JCVI group, create a dataset of:
## the brands to compare and tthe number of comparisons for each brand
brand_comparisons <- second_vax_period_dates %>%
  distinct(brand, n_comparisons) %>%
  group_by(brand) %>%
  mutate(n = n()) %>%
  ungroup()

# check that this is as expected
condition <- 
  n_distinct(brand_comparisons$brand) %in% c(1,2) && 
  max(brand_comparisons$n) == 1
stopifnot(
  "There must be one or two brands, and one value of K per brand." = condition)

# derive comparison arms for k comparisons ----
# define start_fu_date & end_fu_date for each comparison
comparison_arms <- function(
  b, # vaccine brand
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
  
  # vaccinated arm for brand b comparison k
  data_vax <- data_eligible_c %>%
    filter(brand %in% b) %>%
    # start date for vax arm depends on second vax date
    mutate(start_fu_date = covid_vax_2_date + days(d)) %>%
    # no third dose before start_fu_date
    filter(no_evidence_of(covid_vax_3_date, start_fu_date)) %>%
    mutate(arm = "vax") %>%
    droplevels()
  
  # unvaccinated arm for brand b comparison k
  data_unvax <- data_eligible_d %>%
    filter(brand %in% b) %>%
    # start_fu_date for unvax arm depends on elig_date, region and brand
    mutate(start_fu_date = start_of_period + days(d)) %>%
    filter(
      # no first dose before start_fu_date
      no_evidence_of(covid_vax_1_date, start_fu_date),
      # only keep 50% of unvax individuals, depending on if k odd or even
      split %in% split_string
    ) %>%
    mutate(arm = "unvax") %>%
    select(-split) %>%
    droplevels()
  
  # bind datasets from both arms and apply exclusions
  out <- bind_rows(data_vax, data_unvax) %>%
    left_join(data_covs %>%
                select(patient_id, 
                       all_of(exclude_if_evidence_of)),
              by = "patient_id") %>%
    # exclude if evidence of xxx before start_fu_date
    filter_at(all_of(exclude_if_evidence_of),
              all_vars(no_evidence_of(., start_fu_date))) %>%
    select(patient_id, jcvi_group, elig_date, region_0, brand, arm, start_fu_date) 
  
  # elig_date:brand:region-specific end_fu_date for unvax arm
  end_fu_dates <- out %>%
    group_by(jcvi_group, elig_date, brand, region_0) %>%
    # each unvaxxed individual followed up for 2*28 days
    summarise(end_fu_date = min(start_fu_date) + 2*days(28), 
              .groups = "keep") %>%
    ungroup() %>%
    mutate(arm = "unvax")
  
  # join and derive individual-specific end_fu_date for vax arm
  out <- out %>%
    left_join(end_fu_dates, 
              by = c("jcvi_group", "elig_date", "region_0", "brand", "arm")) %>%
    mutate(across(end_fu_date,
                  ~ if_else(arm %in% "vax", 
                            # each individual in vax arm followed up for 28 days
                            start_fu_date + days(28), 
                            .x))) %>% 
    mutate(comparison = k)
  
  # check follow up times are correct in the two arms
  check_fu_time <- out %>%
    mutate(fu_time = as.numeric(end_fu_date - start_fu_date)) %>%
    group_by(arm) %>%
    summarise(min_fu = min(fu_time), max_fu = max(fu_time)) %>%
    ungroup() %>%
    mutate(check = if_else(
      arm == "unvax",
      min_fu == 56 & min_fu == max_fu,
      min_fu == 28 & min_fu == max_fu
    ))
  stopifnot(
    "Follow-up times must be 28 days in vax arm and 56 days in unvax arm" = all(check_fu_time$check))
  
  return(out)
  
}

data_comparison_arms <- bind_rows(lapply(
  as.character(brand_comparisons$brand),
  function(x)
    bind_rows(
      lapply(
        1:brand_comparisons$n_comparisons[brand_comparisons$brand == x],
        function(y)
          comparison_arms(b = x, k = y))))) %>%
  mutate(across(arm, factor, levels = c("vax", "unvax"))) %>%
  mutate(across(comparison, factor))

# check that no overlap between individuals in:
## vax and unvax arms
## unvax arm in odd and even comparisons
check_overlap <- data_comparison_arms %>%
  select(patient_id, arm, comparison) %>%
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
  distinct(patient_id, name, value) %>%
  pivot_wider(
    names_from = name, values_from = value) %>%
  mutate(
    vax_unvax = is.na(vax + unvax_even) & is.na(unvax_odd + vax),
    odd_even = is.na(unvax_odd + unvax_even))

stopifnot("Overlap between vax and unvax arms" = all(check_overlap$vax_unvax))
stopifnot("Overlap between odd and even splits" = all(check_overlap$odd_even))
         
################################################################################
# read long datasets from recurring variables ----

# shielded (index is time_zero)
data_shielded <- data_comparison_arms %>%
  distinct(patient_id, comparison, brand, start_fu_date) %>%
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
  filter(date <= start_fu_date) %>%
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
  distinct(patient_id, brand, comparison, start_fu_date) %>%
  inner_join(readr::read_rds(
    here::here("output", "data", "data_long_bmi_dates.rds")) %>%
      select(patient_id, date, bmi),
            by = "patient_id") %>%
  filter(date <= start_fu_date ) %>%
  group_by(patient_id, brand, comparison) %>%
  arrange(desc(date), .by_group = TRUE) %>%
  # keeps just the most recent before start_fu_date
  distinct(patient_id, comparison, .keep_all = TRUE) %>%
  ungroup() %>%
  select(patient_id, brand, comparison, bmi)

# hospital admissions
# TODO (maybe)

################################################################################

# add clinical and demographic covariates to data_comparison_arms ----

strata_vars <- c("region", "elig_date", "comparison")
demographic_vars <- c("age", "sex", "ethnicity", "imd_0")
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
    data_covs %>%
      select(patient_id, 
             dob, 
             all_of(demographic_vars[demographic_vars %in% names(.)]),
             all_of(ever_vars),
             all_of(clinical_vars[clinical_vars %in% names(.)]),
             all_of(end_vars)), 
    by = "patient_id") %>%
  mutate(imd = factor(imd_0)) %>%
  left_join(
    data_shielded, 
    by = c("patient_id", "brand", "comparison")) %>%
  mutate(across(
    shielded, 
    ~ if_else(is.na(.x), FALSE, .x))) %>%
  left_join(
    data_bmi, 
    by = c("patient_id", "brand", "comparison")) %>%
  mutate(across(bmi, as.character)) %>%
  mutate(across(bmi,
                ~ factor(
                  if_else(is.na(.x), "Not obese", .x),
                  levels = levels(data_bmi$bmi)))) %>%
  mutate(across(
    all_of(ever_vars), 
     ~ case_when(
       is.na(.x) ~ FALSE,
       .x <= start_fu_date ~ TRUE,
       TRUE ~ FALSE))) %>%
  rename_with(
    .fn = ~ str_remove(.x, "_date"),
    .cols = all_of(ever_vars)
    ) %>%
  mutate(
    
    age = as.numeric(start_fu_date - dob)/365.25,
    
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
   patient_id, elig_date, region = region_0, ethnicity, brand, arm, 
   start_fu_date, end_fu_date, comparison,
   dob, 
   unname(unlist(model_varlist))
  ) %>%
  droplevels()

readr::write_rds(
  data_comparisons,
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds"),
  compress = "gz")
