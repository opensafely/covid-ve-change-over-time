library(tidyverse)
library(survival)
library(broom)
library(broom.helpers)

################################################################################
# define formulas

formula_unadj <- Surv(tstart, tstop, status) ~ comparison:arm + strata(region) + strata(comparison)
formula_demog <- . ~ . + poly(age, degree=2, raw=TRUE) + sex + imd + ethnicity
formula_clinical <- . ~ . +
  
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
  formula_unadj,
  formula_demog,
  formula_clinical,
)


model_names = c(
  "Unadjusted" = "0",
  "Adjusting demographics" = "1",
  "Adjusting for demographics + clinical" = "2"
)

formula0 <- formula_unadj
formula1 <- formula_unadj %>% update(formula_demog)
formula2 <- formula_unadj %>% update(formula_demog) %>% update(formula_clinical)


################################################################################
# define model

opt_control <- coxph.control(iter.max = 30)

cox_model <- function(
  number, 
  formula_cox, 
  filename_prefix
  ) {
  
  coxmod <- coxph(
    formula = formula_cox,
    data = data_cox,
    robust = TRUE,
    id = patient_id,
    na.action = "na.fail",
    control = opt_control
  )
  
  # print(warnings())
  # logoutput(
  #   glue("model{number} data size = ", coxmod$n),
  #   glue("model{number} memory usage = ", format(object.size(coxmod), units="GB", standard="SI", digits=3L)),
  #   glue("convergence status: ", coxmod$info[["convergence"]])
  # )
  
  tidy <-
    broom.helpers::tidy_plus_plus(
      coxmod,
      exponentiate = FALSE
    ) %>%
    add_column(
      model = number,
      .before=1
    )
  
  glance <-
    broom::glance(coxmod) %>%
    add_column(
      model = number,
      convergence = coxmod$info[["convergence"]],
      ram = format(object.size(coxmod), units="GB", standard="SI", digits=3L),
      .before = 1
    )
  
  coxmod$data <- NULL
  
  readr::write_rds(
    coxmod,
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{filename_prefix}_model{number}.rds")), 
    compress="gz")
  
  lst(glance, tidy)
}
