################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  group <- "02"
  
} else{
  group <- args[[1]]
}

################################################################################
second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds")) %>%
  filter(jcvi_group %in% group, include) %>%
  distinct(brand, n_comparisons)


outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)

################################################################################
# function for printing glance
print_table <- function(b, outcome) {
  
  model_glance <- readr::read_csv(here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_glance.csv")))
  
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

for (b in as.character(unique(second_vax_period_dates$brand))) {
  capture.output(
        map(outcomes, function(x) try(print_table(b, outcome = x))),
        file = here::here("output", glue("jcvi_group_{group}"), "tables", glue("{b}_modelcox_glance.txt")),
        append=FALSE
      )
}

################################################################################
# function for printing coefficients from all three models
print_coefficients <- function(b, outcome) {
  
  data <- readr::read_rds(
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_summary.rds")))
               
  data %>%
    filter(outcome %in% outcome, !str_detect(term, "unvax")) %>%
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

for (b in as.character(unique(second_vax_period_dates$brand))) {
  capture.output(
    map(outcomes, function(x) try(print_coefficients(b, outcome = x))),
    file = here::here("output", glue("jcvi_group_{group}"), "tables", glue("{b}_modelcox_coefficients.txt")),
    append=FALSE
  )
}


