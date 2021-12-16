################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobs <- FALSE
  group <- "02"
  outcome <- "postest"
  
} else{
  removeobs <- TRUE
  group <- args[[1]]
  outcome <- args[[2]]
}

################################################################################
# read data

data_comparisons <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds"))

source(here::here("analysis", "lib", "model_cox.R"))

model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

################################################################################
# apply model


for (b in unique(data_comparisons$brand)) {
 
  # read data_tte
  data_tte <- readr::read_rds(
    here::here("output", glue("jcvi_group_{group}"), "data", glue("data_tte_{b}_{outcome}.rds")))

  data_cox <- data_tte %>%
    left_join(data_comparisons %>%
                filter(brand %in% b) %>%
                select(patient_id, comparison, elig_date, region, 
                       unname(unlist(model_varlist))),
              by = c("patient_id", "comparison")) 
  
  # apply model  
  summary0 <- cox_model(
    number = 0, 
    formula_cox = formula0, 
    filename_prefix = glue("{b}_{outcome}"))
  summary1 <- cox_model(
    number = 1, 
    formula_cox = formula1, 
    filename_prefix = glue("{b}_{outcome}"))
  summary2 <- cox_model(
    number = 2, 
    formula_cox = formula2, 
    filename_prefix = glue("{b}_{outcome}"))
  
  
  # combine results
  model_glance <-
    bind_rows(summary0$glance, summary1$glance, summary2$glance) %>%
    mutate(
      model_name = fct_recode(as.character(model), !!!model_names),
      outcome = outcome
    )
  
  readr::write_csv(
    model_glance,
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_glance.csv"))) 
  
  model_tidy <- bind_rows(summary0$tidy, summary1$tidy, summary2$tidy) %>%
    mutate(
      model_name = fct_recode(as.character(model), !!!model_names),
      outcome = outcome
    )
  
  readr::write_rds(
    model_tidy,
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_tidy.rds")), 
    compress="gz")
  readr::write_csv(
    model_tidy,
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_tidy.csv")))
  

}

