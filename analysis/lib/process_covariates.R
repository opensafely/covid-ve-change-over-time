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

data_tests <- arrow::read_feather(
  file = here::here("output", "input_tests.feather")) %>%
  transmute(
    patient_id, 
    covid_test_n = covid_test_pre_elig_n + covid_test_post_elig_n
    )

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
  demographic_vars <- c("age_band", "sex", "ethnicity", "imd")
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
    "bmi",
    "covid_test_n"
  )
  
  out <- .data %>%
    # join and process covariates
    left_join(
      data_processed %>%
        select(patient_id, subgroup, 
               all_of(censor_vars),
               all_of(demographic_vars[demographic_vars %in% names(.)]),
               all_of(ever_vars),
               all_of(clinical_vars[clinical_vars %in% names(.)]),
               all_of(str_c(outcomes, "_date"))), 
      by = "patient_id") %>%
    mutate(imd = factor(imd)) %>%
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
      patient_id, jcvi_group, elig_date, region, ethnicity, arm, subgroup,
      start_fu_date, end_fu_date, comparison,
      unname(unlist(model_varlist)),
      all_of(str_c(outcomes, "_date")),
      all_of(censor_vars)
    ) %>%
    droplevels()
  
  return(out)
  
}
