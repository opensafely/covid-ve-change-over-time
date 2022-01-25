################################################################################

# This script:
# applies the cox model for a given comparison and outcome


################################################################################
library(tidyverse)
library(fastDummies)
library(glue)
library(survival)
library(broom)
library(broom.helpers)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "BNT162b2"
  subgroup_label <- 1
  outcome <- "postest"
  
} else{
  comparison <- args[[1]]
  subgroup_label <- args[[2]]
  outcome <- args[[3]]
}

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroups <- c(subgroups, "all")
subgroup <- subgroups[subgroup_label]

################################################################################
# read data

model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

data_cox <- readr::read_rds(
  here::here("output", "models_cox", "data", glue("data_cox_{comparison}_{subgroup_label}_{outcome}.rds")))

################################################################################
# define formulas

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
  
  # opt_control <- coxph.control(iter.max = 30)
  if (filename_prefix %in% c("ChAdOx_3_noncoviddeath", "ChAdOx_3_death")) {
    # because adjusted model failing despite no issues with low event counts
    opt_control <- coxph.control(iter.max = 50)
  } else {
    opt_control <- coxph.control(iter.max = 30)
  }
  
  cat(glue("...... fitting model {number} ......"), "\n")
  cat(glue("{model_name}"), "\n")
  timetofit <- system.time((
    coxmod <- coxph(
      formula = formula_cox,
      data = data_cox_dummy,
      robust = TRUE,
      id = patient_id,
      na.action = "na.fail",
      control = opt_control)
  ))
  
  # process model (as having issues with broom tidy)
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
    here::here("output", "models_cox", "data", glue("model{number}_{filename_prefix}.rds")), 
    compress="gz")
  
  lst(glance = glance, summary = coxmod_summary)
}

################################################################################

model_output <- list()
model_output[[1]] <- try(cox_model(
  number = 0, 
  formula = formula_cox_0,
  filename_prefix = glue("{comparison}_{subgroup_label}_{outcome}")))
# discard demographic only adjusted model for now
# model_output[[2]] <- try(cox_model(
#   number = 1,
#   formula = formula_cox_1,
#   filename_prefix = glue("{subgroup}_{comparison}_{outcome}")))
model_output[[3]] <- try(cox_model(
  number = 2, 
  formula = formula_cox_2,
  filename_prefix = glue("{comparison}_{subgroup_label}_{outcome}")))


# check for errors 
check_errors <- sapply(model_output, function(x) any(class(x) %in% "try-error"))

if (!all(check_errors)) {
  
  model_summary <- bind_rows(
    lapply(
      # only bind tibbles (to avoid errors in case some models did not converge)
      seq_along(model_output)[sapply(model_output, function(x) is_tibble(x[[1]]))],
      # select summary
      function(x) model_output[[x]]$summary
    )) %>%
    mutate(outcome = outcome)
  readr::write_rds(
    model_summary,
    here::here("output", "models_cox", "data", glue("modelcox_summary_{comparison}_{subgroup_label}_{outcome}.rds"))) 
  
  ### postprocessing using broom (may be unreliable)
  # combine results
  model_glance <- bind_rows(
    lapply(
      # only bind tibbles (to avoid errors in case some models did not converge)
      seq_along(model_output)[sapply(model_output, function(x) is_tibble(x[[1]]))],
      # select glance
      function(x) model_output[[x]]$glance
    )) %>%
    mutate(outcome = outcome)
  
  readr::write_rds(
    model_glance,
    here::here("output", "models_cox", "data", glue("modelcox_glance_{comparison}_{subgroup_label}_{outcome}.rds"))) 
  
}
