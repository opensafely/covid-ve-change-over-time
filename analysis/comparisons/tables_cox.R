################################################################################

# This script:
# - creates and saves tables of metadata and coefficients from the cox models

################################################################################
library(tidyverse)
library(glue)
## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "BNT162b2"
  subgroup_label <- 2
  
} else{
  comparison <- args[[1]]
  subgroup_label <- args[[2]]
}

################################################################################
fs::dir_create(here::here("output", "models_cox",  "tables"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "second_vax_period", "data", "second_vax_period_dates.rds")) 

outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)

subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds")
)

################################################################################
# function for printing glance
print_table <- function(outcome) {
  
  model_glance <- readr::read_rds(
    here::here("output", "models_cox", "data", glue("modelcox_glance_{comparison}_{subgroup_label}_{outcome}.rds")))
  
  model_glance %>%
    select(-outcome) %>%
    mutate(across(where(is.numeric), signif, digits = 6)) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(
      cols = -model
    ) %>% 
    pivot_wider(
      names_from = "model",
      values_from = value
    ) %>%
    kableExtra::kable("pipe", caption = outcome)
  
}


capture.output(
  map(outcomes, function(x) try(print_table(outcome = x))),
  file = here::here("output", "models_cox", "tables", glue("modelcox_glance_{comparison}_{subgroup_label}.txt")),
  append=FALSE
)


################################################################################
# function for printing coefficients from all three models
print_coefficients <- function(outcome) {
  
  data <- readr::read_rds(
    here::here("output", "models_cox", "data", glue("modelcox_summary_{comparison}_{subgroup_label}_{outcome}.rds")))
               
  data %>%
    filter(outcome %in% outcome) %>%
    mutate(across(c(estimate, lower, upper), round, 2)) %>%
    mutate(value = str_c(estimate, " (", lower, "-", upper, ")")) %>%
    select(term, value, model) %>% 
    mutate(across(model, 
                  factor,
                  levels = as.character(0:2),
                  labels = c("unadjusted",
                             "demographic",
                             "demographic + clinical"))) %>%
    pivot_wider(names_from = "model", values_from = "value") %>%
    kableExtra::kable("pipe", caption = outcome)
  
}

capture.output(
  map(outcomes, function(x) try(print_coefficients(outcome = x))),
  file = here::here("output", "models_cox", "tables", glue("modelcox_coefficients_{comparison}_{subgroup_label}.txt")),
  append=FALSE
)


