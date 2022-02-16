################################################################################
# process comparisons data
library(tidyverse)
library(lubridate)
library(glue)

################################################################################
## create folders for outputs
fs::dir_create(here::here("output", "data"))

################################################################################
## import study_parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))
K <- study_parameters$max_comparisons

# import variable names
model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

data_eligible_e <- readr::read_csv(
  here::here("output", "data", "data_eligible_e.csv"))

# read data for ever covariates
data_ever <- arrow::read_feather(
  file = here::here("output", "input_ever.feather")) 
# read data for k covariates
data_k <- bind_rows(lapply(
  1:2,
  function(k)
  arrow::read_feather(
    file = here::here("output", glue("input_{k}.feather"))) %>%
    mutate(k=k)
))

################################################################################
# process covaraiets data
data_covariates <- data_k %>% 
  left_join(data_ever %>% select(-start_1_date), 
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
  # define clinical covariates
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
    ),
    # chronic respiratory disease other than asthma ever
    other_respiratory = resp_date <= start_k_date,
    chronic_respiratory_disease = asthma | other_respiratory,
    # chronic neurological disease including significant learning disorder
    chronic_neuro_inc_ld = cns_date <= start_k_date,
    # wider learning disorder
    ld_inc_ds_and_cp = learndis_date <= start_k_date,
    # diabetes ever
    diabetes = diab_date <= start_k_date,
    # severe mental illness ever
    sev_ment = sev_mental_date <= start_k_date,
    # chronic heart disease ever
    chronic_heart_disease = chd_date <= start_k_date,
    # chronic liver disease ever
    chronic_liver_disease = cld_date <= start_k_date,
    # permanent immunosupression
    permanant_immunosuppression = immdx_date <= start_k_date,
    # current immunosuppression medication
    immunosuppression_meds = immrx,
    # asplenia or dysfunction of the spleen ever
    asplenia_ever = spln_date <= start_k_date,
    # chronic kidney disease stages 3-5
    ckd = ckd_group,
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
      ),
    
    pregnancy = preg_group
  ),
  
  any_immunosuppression = (
    permanant_immunosuppression | 
      asplenia | 
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
  select(patient_id, start_k_date, end_k_date, 
         all_of(unlist(unname(model_varlist))))
  

readr::write_rds(
  data_covariates,
  here::here("output", "data", "data_covariates.rds"),
  compress = "gz"
)

