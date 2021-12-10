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

