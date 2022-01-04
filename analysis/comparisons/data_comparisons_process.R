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
  K = study_parameters$max_comparisons 
) {
  
  comparison_k <- function(k) {
    
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
    out_k <- data %>%
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
    "efi",
    "bmi"
  )
  end_vars <- c(
    "coviddeath_date",
    "death_date",
    "dereg_date"
  )
  
  # define breaks and labels for age_band
  age_breaks_lower <- c(16, seq(30,90,10))
  age_breaks_upper <- as.character(lead(age_breaks_lower) - 1)
  age_breaks_upper[-length(age_breaks_upper)] <- str_c("-", age_breaks_upper[-length(age_breaks_upper)])
  age_breaks_upper[length(age_breaks_upper)] <- "+"
  age_labels <- str_c(age_breaks_lower, age_breaks_upper)
  
  out <- .data %>%
    mutate(arm = str_c(brand, arm, sep = "_")) %>%
    select(-brand) %>%
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
      
      # change age band labels from "(lower,upper]" to "lower-upper"
      age_band = str_remove(age_band, "\\("),
      age_band = str_remove(age_band, "\\["),
      age_band = str_remove(age_band, "\\]"),
      
      age_band = str_replace(age_band, ",Inf", "+"),
      age_band = str_replace(age_band, ",", "-"),
      age_band = factor(age_band),
      
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
      patient_id, jcvi_group, elig_date, region = region_0, ethnicity, arm, 
      start_fu_date, end_fu_date, comparison,
      dob, 
      unname(unlist(model_varlist))
    ) %>%
    droplevels()
  
  return(out)
  
}

################################################################################
# need to save an unvax dataset corresponding to each brand - even though
# these are the same individuals, the dates on which covariates are defined 
# depend on the brand in the vax arm
for (b in c("BNT162b2", "ChAdOx")) {
  for (a in c("vax", "unvax")) {
    data_comparisons <- comparison_arms(b=b, a=a) %>%
      process_covariates()
    
    readr::write_rds(
      data_comparisons,
      here::here("output", "data", glue("data_comparisons_{b}_{a}.rds")),
      compress = "gz"
    )
    
  }
}

