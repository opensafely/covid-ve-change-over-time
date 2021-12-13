


################################################################################

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  group <- "02"
  outcome <- "postest"
  
} else{
  group <- args[[1]]
  outcome <- args[[2]]
}

# outcomes <- c("postest", "covidadmitted", "coviddeath", "death")
censor <- c("noncoviddeath", "dereg")

data_tte <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"),  "data", "data_tte.rds"))


data_tte <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"),  "data", "data_tte.rds")) %>%
  select(patient_id, elig_date, region, brand, arm, time_zero, end_fu, comparison, 
         starts_with(str_c(c(outcome, censor)))) %>%
  group_by(patient_id, brand, comparison) %>%
  # censor on noncoviddeath or dereg date
  mutate(
    censor_tte = min(noncoviddeath_tte, dereg_tte),
    censor_reason = case_when(
      noncoviddeath_tte == end_fu & dereg_tte == end_fu ~ "nocensor",
      noncoviddeath_tte < dereg_tte ~ "noncoviddeath",
      dereg_tte < noncoviddeath_tte ~ "dereg_date",
      TRUE ~ NA_character_
    )) %>%
  ungroup() %>%
  select(-starts_with(censor))

apply_censoring <- function(.data) {
  
  .data
  
  .data %>%
    mutate(across(glue("{outcome}_status"),
                  if_else(
                    !! glue("{outcome}_tte") < dereg_status
                  ))
    mutate(across(glue("{outcome}_tte"),
                  min()))
    mutate(
      
    )
}


formula_unadj <- Surv(start, end, status) ~ k:arm + strata(elig_date) + strata(region) + strata(k)
formula_clin <- . ~ . + clin


library(survival)
mod_test <- coxph(formula_unadj, data = data_tte)
sum_test <- summary(mod_test)
sum_test$
  


# formulas ----

formula_demog <- . ~ . + poly(age, degree=2, raw=TRUE) + sex + imd + ethnicity
formula_comorbs <- . ~ . +
  
  bmi +
  heart_failure +
  other_heart_disease +
  
  dialysis +
  diabetes +
  chronic_liver_disease +
  
  current_copd +
  #cystic_fibrosis +
  other_respiratory +
  
  lung_cancer +
  haematological_cancer +
  cancer_excl_lung_and_haem +
  
  any_immunosuppression +
  
  dementia +
  other_neuro_conditions +
  
  ld_inc_ds_and_cp +
  psychosis_schiz_bipolar +
  
  multimorb +
  
  shielded +
  
  flu_vaccine +
  
  efi


list_formula <- lst(
  formula_demog,
  formula_comorbs,
)

