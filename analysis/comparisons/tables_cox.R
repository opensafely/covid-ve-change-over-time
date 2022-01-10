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
  subgroup <- "03-10"
  comparison <- "BNT162b2"
  
} else{
  subgroup <- args[[1]]
  comparison <- args[[2]]
}

################################################################################
fs::dir_create(here::here("output", "tables"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds")) 

outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)

################################################################################
# function for printing glance
print_table <- function(outcome) {
  
  model_glance <- readr::read_csv(here::here("output", "models", glue("modelcox_glance_{subgroup}_{comparison}_{outcome}.csv")))
  
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
  file = here::here("output", "tables", glue("modelcox_glance_{subgroup}_{comparison}.txt")),
  append=FALSE
)


################################################################################
# function for printing coefficients from all three models
print_coefficients <- function(outcome) {
  
  data <- readr::read_rds(
    here::here("output", "models", glue("modelcox_summary_{subgroup}_{comparison}_{outcome}.rds")))
               
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
  file = here::here("output", "tables", glue("modelcox_coefficients_{subgroup}_{comparison}.txt")),
  append=FALSE
)


