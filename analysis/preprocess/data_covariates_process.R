################################################################################
# process comparisons data
library(tidyverse)
library(lubridate)
library(glue)

################################################################################
## import study_parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))
K <- study_parameters$max_comparisons

# import variable names
model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

# individuals eligible based on box c, d & e criteria 
# arm and split info
data_arm <- bind_rows(
  readr::read_rds(
    here::here("output", "data", "data_eligible_e_vax.rds")) %>%
    rename(arm=brand),
  readr::read_rds(
    here::here("output", "data", "data_eligible_e_unvax.rds")) %>%
    mutate(arm = "unvax")
)  %>%
  select(patient_id, arm, split)

# sex from data_processed
data_sex <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) %>%
  select(patient_id, sex)

# read data for ever covariates
data_ever <- arrow::read_feather(
  file = here::here("output", "input_ever.feather")) 
# read data for k covariates
data_k <- bind_rows(lapply(
  1:K,
  function(k)
  arrow::read_feather(
    file = here::here("output", glue("input_{k}.feather"))) %>%
    mutate(k=k)
))

ever_before <- function(.data, name, var) {
  .data %>%
    mutate(!! sym(name) := if_else(
      !is.na(!! sym(var)) & (!! sym(var) <= start_k_date),
      TRUE,
      FALSE
    ))
}

################################################################################
# process covariates data
data_covariates <- data_arm %>% 
  # join ever covariates
  left_join(data_ever %>% select(-start_1_date), 
            by = "patient_id") %>%
  # join period-updating covariates
  left_join(data_k,
            by = "patient_id") %>%
  # join sex
  left_join(data_sex,
            by = "patient_id") %>%
  # clean BMI data
  mutate(across(bmi_stage,
                ~ case_when(
                  is.na(.x) | .x %in% "Decreased body mass index"
                  ~ NA_character_,
                  .x %in% c("Body mass index 30+ - obesity",
                            "Obese class I",
                            "Obese class I (body mass index 30.0 - 34.9)")
                  ~ "Obese I (30-34.9)",
                  .x %in% c("Obese class II",
                            "Obese class II (body mass index 35.0 - 39.9)")
                  ~ "Obese II (35-39.9)",
                  .x %in% c("Body mass index 40+ - severely obese",
                            "Obese class III",
                            "Obese class III (body mass index equal to or greater than 40.0)")
                  ~ "Obese III (40+)",
                  TRUE
                  ~ "Not obese"
                ))) %>%
  mutate(across(bmi_stage_date, 
                ~ if_else(
                  is.na(bmi_stage),
                  as.POSIXct(NA_character_),
                  .x))) %>%
  mutate(across(bmi,
                ~ case_when(
                  .x < 10 | bmi >= 100 
                  ~ NA_character_,
                  .x < 30 
                  ~ "Not obese", 
                  .x >= 30 & .x < 35 
                  ~ "Obese I (30-34.9)",
                  .x >= 35 & .x < 40 
                  ~ "Obese II (35-39.9)",
                  .x >= 40 
                  ~ "Obese III (40+)",
                  TRUE ~ NA_character_))) %>%
  mutate(across(bmi_date_measured, 
                ~ if_else(
                  is.na(bmi),
                  as.POSIXct(NA_character_),
                  .x))) %>%
  # clean asthma data
  mutate(
    # clinically extremely vulnerable in period k
    cev = cev_group,
    # poorly controlled asthma in period k
    asthma = case_when(
      astadm ~ TRUE,
      !is.na(astdx_date) & 
        astdx_date <= start_k_date & 
        astrxm1 &
        astrxm2 & 
        astrxm3 ~ TRUE,
      TRUE ~ FALSE
    )) %>%
  # clean test history data
  mutate(across(test_hist_n,
                ~ factor(case_when(
                  is.na(.x) ~ NA_character_,
                  .x < 1 ~ "0",
                  .x < 2 ~ "1",
                  .x < 3 ~ "2",
                  TRUE ~ "3+"
                )))) %>%
  # clean "ever" variables
  # chronic respiratory disease other than asthma ever
  ever_before(
    name = "other_respiratory",
    var = "resp_date"
  ) %>%
  # chronic neurological disease including significant learning disorder
  ever_before(
    name = "chronic_neuro_inc_ld",
    var = "cns_date"
  ) %>%
  # wider learning disorder
  ever_before(
    name = "ld_inc_ds_and_cp",
    var = "learndis_date"
  ) %>%
  # diabetes ever
  ever_before(
    name = "diabetes",
    var = "diab_date"
  ) %>%
  # severe mental illness ever
  ever_before(
    name = "sev_ment",
    var = "sev_mental_date"
  ) %>%
  # chronic heart disease ever
  ever_before(
    name = "chronic_heart_disease",
    var = "chd_date"
  ) %>%
  # chronic liver disease ever
  ever_before(
    name = "chronic_liver_disease",
    var = "cld_date"
  ) %>%
  # permanent immunosupression
  ever_before(
    name = "permanant_immunosuppression",
    var = "immdx_date"
  ) %>%
  # asplenia or dysfunction of the spleen ever
  ever_before(
    name = "asplenia_ever",
    var = "spln_date"
  ) %>%
  # resident in longterm residential home
  ever_before(
    name = "longres",
    var = "longres_date"
  ) %>%
  mutate(
    # current immunosuppression medication
    immunosuppression_meds = immrx,
    # chronic kidney disease stages 3-5
    ckd = ckd_group,
    # any chronic respiratory disease
    chronic_respiratory_disease = asthma | other_respiratory,
    # bmi
    bmi = factor(
      case_when(
        is.na(bmi) & is.na(bmi_stage) ~ "Not obese",
        is.na(bmi) ~ bmi_stage,
        is.na(bmi_stage) ~ bmi,
        bmi_stage_date <= bmi_date_measured ~ bmi,
        TRUE ~ bmi_stage
    ),
    levels = c(
      "Not obese",
      "Obese I (30-34.9)",
      "Obese II (35-39.9)",
      "Obese III (40+)"
      )
  ),
  
  pregnancy = preg_group & (sex == "Female") & (age < 50),
  
  any_immunosuppression = (
    permanant_immunosuppression | 
      asplenia_ever | 
      immunosuppression_meds),
  
  multimorb =
    bmi %in% "Obese III (40+)" +
    chronic_heart_disease  +
    diabetes +
    chronic_liver_disease +
    ckd +
    chronic_respiratory_disease +
    any_immunosuppression +
    chronic_neuro_inc_ld +
    ld_inc_ds_and_cp +
    sev_ment,
  multimorb = cut(
    multimorb, 
    breaks = c(0, 1, 2, Inf),
    labels=c("0", "1", "2+"), 
    right=FALSE)
  
  ) %>%
  mutate(across(contains("_date"), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days"))) %>%
  select(patient_id, start_k_date, end_k_date, k,
         arm, split,
         anytest_date, age,
         all_of(unname(model_varlist$clinical)))
  

readr::write_rds(
  data_covariates,
  here::here("output", "data", "data_covariates.rds"),
  compress = "gz"
)


