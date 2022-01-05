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
  outcome <- "postest"
  
} else{
  comparison <- args[[1]]
  outcome <- args[[2]]
}

################################################################################
# read data

fs::dir_create(here::here("output", "models"))
fs::dir_create(here::here("output", "tables"))

model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

################################################################################

if (comparison %in% c("BNT162b2", "ChAdOx")) {
  
  trt <- str_c(comparison, "_vax")
  
  data <- readr::read_rds(
    here::here("output", "data", glue("data_tte_{comparison}_{outcome}.rds"))) %>%
    left_join(
      bind_rows(
        lapply(
          c("vax", "unvax"),
          function(x)
            readr::read_rds(
              here::here("output", "data", glue("data_comparisons_{comparison}_{x}.rds"))) %>%
            select(patient_id, comparison, jcvi_group, elig_date, region, 
                   unname(unlist(model_varlist)))
        )
      ),
      by = c("patient_id", "comparison"))
    
        
} else if (comparison %in% "both") {
  
  trt <- "BNT162b2_vax"
  
  data <- readr::read_rds(
    here::here("output", "data", glue("data_tte_{comparison}_{outcome}.rds"))) %>%
    left_join(
      bind_rows(
        lapply(
          c("BNT162b2", "ChAdOx"),
          function(x)
            readr::read_rds(
              here::here("output", "data", glue("data_comparisons_{x}_vax.rds"))) %>%
            select(patient_id, comparison, jcvi_group, elig_date, region, 
                   unname(unlist(model_varlist)))
        )
      ),
      by = c("patient_id", "comparison"))
  
}

################################################################################

data_strata <- data %>% 
  mutate(strata_var = factor(str_c(jcvi_group, elig_date, region, sep = ", ")))

################################################################################

events_per_personyears <- function(strata) {
  data_strata %>%
    filter(strata_var %in% strata) %>%
    mutate(days = tstop-tstart) %>%
    group_by(comparison, arm) %>%
    summarise(
      n = n(),
      personyears = sum(days)/365,
      events = sum(status),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    mutate(across(c(n, events, personyears),
                  # round to nearest 10
                  ~ scales::comma(round(.x, -1), accuracy = 1))) %>%
    mutate(value = str_c(events, " / ", personyears, " (",n,")")) %>%
    select(comparison, arm, value) %>%
    pivot_wider(names_from = arm, values_from = value) %>%
    kableExtra::kable(
      "pipe",
      caption = glue("Strata: {strata}; {outcome} events / person-years (n)")
    )
}

################################################################################

capture.output(
  for (s in levels(data_strata$strata_var)) {
    print(events_per_personyears(s))
  },
  file = here::here("output", "tables", glue("{comparison}_{outcome}_incidence.txt")),
  append=FALSE
)

################################################################################

data_cox <- data_strata %>%
  dummy_cols(
    select_columns = c("comparison"),
    remove_selected_columns = TRUE
  ) %>%
  mutate(across(starts_with("comparison"),
                ~ if_else(arm %in% trt,
                          .x, 0L))) %>%
  rename_at(vars(starts_with("comparison")), ~str_c(.x, "_", str_remove(trt, "_vax")))

################################################################################
# define formulas

comparisons <- data_cox %>% select(starts_with("comparison")) %>% names()

formula_unadj <- as.formula(str_c(
  "Surv(tstart, tstop, status, type = \"counting\") ~ ",
  str_c(comparisons, collapse = " + "),
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
  
  shielded +
  
  flu_vaccine #+
  
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
    here::here("output", "models", glue("{filename_prefix}_model{number}.rds")), 
    compress="gz")
  
  lst(glance = glance, summary = coxmod_summary)
}

################################################################################

model_output <- list()
model_output[[1]] <- cox_model(
  number = 0, 
  formula = formula_cox_0,
  filename_prefix = glue("{comparison}_{outcome}"))
model_output[[2]] <- cox_model(
  number = 1, 
  formula = formula_cox_1,
  filename_prefix = glue("{comparison}_{outcome}"))
model_output[[3]] <- cox_model(
  number = 2, 
  formula = formula_cox_2,
  filename_prefix = glue("{comparison}_{outcome}"))


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
  here::here("output", "models", glue("{comparison}_{outcome}_modelcox_summary.rds"))) 

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

readr::write_csv(
  model_glance,
  here::here("output", "models", glue("{comparison}_{outcome}_modelcox_glance.csv"))) 

