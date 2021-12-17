library(tidyverse)
library(survival)
library(broom)
library(broom.helpers)
library(tictoc)

################################################################################
# define model

cox_model <- function(
  number, 
  filename_prefix
  ) {
  
  # so that when cox_model is run in lapply the options are printed after each
  # evaluation and not all at the end
  op_warn <- options("warn")
  on.exit(options(op_warn))
  
  options(warn=1)
  
  number <- as.character(number)
  
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
  
  
  model_names = c(
    "0" = "Unadjusted",
    "1" = "Adjusting demographics",
    "2" = "Adjusting for demographics + clinical"
  )
  
  if (number == "0") {
    formula_cox <- formula_unadj
  } else if (number == "1") {
    formula_cox <- formula_unadj %>% update(formula_demog)
  } else if (number == "2") {
    formula_cox <- formula_unadj %>% update(formula_demog) %>% update(formula_clinical)
  }
  
  model_name <- model_names[number]
  
  opt_control <- coxph.control(iter.max = 30)
  
  cat(glue("...... fitting model {number} ......"), "\n")
  cat(glue("{model_name}"), "\n")
  # tic("time to fit model")
  timetofit <- system.time((
    coxmod <- coxph(
      formula = formula_cox,
      data = data_cox,
      robust = TRUE,
      id = patient_id,
      na.action = "na.fail",
      control = opt_control)
  ))
  # toc()
  
  
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
    ) %>%
    mutate(across(
      model,
      factor, levels = names(model_names), labels = model_names)) %>%
    # add output of system.time
    bind_cols(as_tibble(t(as.matrix(timetofit))))
  
  coxmod$data <- NULL
  
  readr::write_rds(
    coxmod,
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{filename_prefix}_model{number}.rds")), 
    compress="gz")
  
  lst(glance, tidy)
}
