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
        patient_id, jcvi_group, elig_date, region_0, arm, 
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
data_long_shielded_dates <- readr::read_rds(
  here::here("output", "data", "data_long_shielded_dates.rds")) %>%
  select(patient_id, date) %>%
  mutate(name = "shielded")

data_long_nonshielded_dates <- readr::read_rds(
  here::here("output", "data", "data_long_nonshielded_dates.rds")) %>%
  select(patient_id, date) %>%
  mutate(name = "nonshielded")

data_long_bmi_dates <- readr::read_rds(
  here::here("output", "data", "data_long_bmi_dates.rds")) %>%
  select(patient_id, date, bmi)

process_covariates <- function(.data) {
  
  # shielded (index is time_zero)
  data_shielded <- .data %>%
    distinct(patient_id, comparison, start_fu_date) %>%
    inner_join(
      bind_rows(data_long_shielded_dates, data_long_nonshielded_dates),
      by = "patient_id") %>%
    filter(date <= start_fu_date) %>%
    group_by(patient_id, comparison, name) %>%
    summarise(date = max(date), .groups = "keep") %>%
    ungroup() %>%
    pivot_wider(names_from = name, values_from = date) %>%
    mutate(shielded = case_when(
      !is.na(shielded) & 
        ((!is.na(nonshielded) & (shielded > nonshielded)) |
           is.na(nonshielded)) ~ TRUE,
      TRUE ~ FALSE)) %>%
    select(patient_id, comparison, shielded)
  
  
  # bmi (index is time_zero)
  data_bmi <- .data %>%
    distinct(patient_id, comparison, start_fu_date) %>%
    inner_join(data_long_bmi_dates, by = "patient_id") %>%
    filter(date <= start_fu_date ) %>%
    group_by(patient_id, comparison) %>%
    arrange(desc(date), .by_group = TRUE) %>%
    # keeps just the most recent before start_fu_date
    distinct(patient_id, comparison, .keep_all = TRUE) %>%
    ungroup() %>%
    select(patient_id, comparison, bmi)
  
  # hospital admissions
  # TODO (maybe)
  
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
    # "efi",
    "bmi"
  )
  
  # define breaks and labels for age_band
  age_breaks_lower <- c(16, seq(20,95,5))
  age_breaks_upper <- as.character(lead(age_breaks_lower) - 1)
  age_breaks_upper[-length(age_breaks_upper)] <- str_c("-", age_breaks_upper[-length(age_breaks_upper)])
  age_breaks_upper[length(age_breaks_upper)] <- "+"
  age_labels <- str_c(age_breaks_lower, age_breaks_upper)
  
  out <- .data %>%
    # join and process covariates
    left_join(
      data_processed %>%
        select(patient_id, 
               dob, subgroup, 
               all_of(censor_vars),
               all_of(demographic_vars[demographic_vars %in% names(.)]),
               all_of(ever_vars),
               all_of(clinical_vars[clinical_vars %in% names(.)]),
               all_of(str_c(outcomes, "_date"))), 
      by = "patient_id") %>%
    mutate(imd = factor(imd_0)) %>%
    left_join(
      data_shielded, 
      by = c("patient_id", "comparison")) %>%
    mutate(across(
      shielded, 
      ~ if_else(is.na(.x), FALSE, .x))) %>%
    left_join(
      data_bmi, 
      by = c("patient_id", "comparison")) %>%
    mutate(across(bmi, as.character)) %>%
    mutate(across(bmi,
                  ~ factor(
                    if_else(is.na(.x), "Not obese", .x),
                    levels = levels(data_bmi$bmi)))) %>%
    # only keep outcome or censor data in the comparison period
    mutate(across(
      all_of(c(censor_vars, str_c(outcomes, "_date"))), 
      ~ case_when(
        is.na(.x) ~ as.Date(NA_character_),
        .x <= start_fu_date ~ as.Date(NA_character_),
        .x <= end_fu_date ~ .x,
        TRUE ~  as.Date(NA_character_)))) %>%
    # process ever covariates
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
      
      # 10-year age bands (apart from lowest, which is 14-year)
      age = floor(as.numeric(start_fu_date - dob)/365.25),
      age_band = cut(
        age,
        breaks = c(age_breaks_lower, Inf),
        include.lowest = TRUE,
        labels = age_labels
      ),
      
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
      
      flu_vaccine = flu_vaccine == 1
      
    ) %>%
    select(
      patient_id, jcvi_group, elig_date, region = region_0, ethnicity, arm, subgroup,
      start_fu_date, end_fu_date, comparison, dob, 
      unname(unlist(model_varlist)),
      all_of(str_c(outcomes, "_date")),
      all_of(censor_vars)
    ) %>%
    droplevels()

  return(out)
  
}

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

