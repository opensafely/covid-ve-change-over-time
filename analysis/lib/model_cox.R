library(tidyverse)
library(survival)
library(broom)
library(broom.helpers)


# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

################################################################################
# define formulas

formula_unadj <- as.formula(str_c(
  "Surv(tstart, tstop, status, type = \"counting\") ~ ",
  str_c(str_c("comprison_", 1:study_parameters$max_comparisons), collapse = " + "),
  " + strata(strata_var)"))
formula_demog <- . ~ . + age_band + sex + imd + ethnicity
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
  
  shielded #+

# flu_vaccine +

# efi


model_names = c(
  "0" = "Unadjusted",
  "1" = "Adjusting demographics",
  "2" = "Adjusting for demographics + clinical"
)

formula_cox_0 <- formula_unadj
formula_cox_1 <- formula_unadj %>% update(formula_demog)
formula_cox_2 <- formula_unadj %>% update(formula_demog) %>% update(formula_clinical)

################################################################################
# define model

cox_model <- function(
  number, 
  formula_cox,
  filename_prefix
  ) {
  
  # so that when cox_model is run in lapply the options are printed after each
  # evaluation and not all at the end
  op_warn <- options("warn")
  on.exit(options(op_warn))
  
  options(warn=1)
  
  number <- as.character(number)
  
  model_name <- model_names[number]
  
  opt_control <- coxph.control(iter.max = 30)
  
  cat(glue("...... fitting model {number} ......"), "\n")
  cat(glue("{model_name}"), "\n")
  timetofit <- system.time((
    coxmod <- coxph(
      formula = formula_cox,
      data = data_cox,
      robust = TRUE,
      id = patient_id,
      na.action = "na.fail",
      control = opt_control)
  ))
  
  # process model (as having issues with broom)
  coxmod_summary <- as_tibble(
    summary(coxmod)$conf.int, 
    rownames = "term"
  ) %>% 
    select(term,
           estimate = `exp(coef)`, 
           lower = `lower .95`, 
           upper = `upper .95`) %>%
    mutate(model = number) 
  
  
  glance <-
    broom::glance(coxmod) %>%
    add_column(
      model = number,
      convergence = coxmod$info[["convergence"]],
      ram = format(object.size(coxmod), units="GB", standard="SI", digits=3L),
      .before = 1
    ) %>%
    mutate(across(
      model,
      factor, levels = names(model_names), labels = model_names)) %>%
    # add output of system.time
    bind_cols(as_tibble(t(as.matrix(timetofit))))
  
  coxmod$data <- NULL
  
  readr::write_rds(
    coxmod,
    here::here("output", "models", glue("{filename_prefix}_model{number}.rds")), 
    compress="gz")
  
  lst(glance = glance, summary = coxmod_summary)
}
