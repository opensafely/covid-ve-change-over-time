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
  subgroup_label <- 1
  
} else{
  comparison <- args[[1]]
  subgroup_label <- args[[2]]
}

################################################################################
fs::dir_create(here::here("output", "models_cox",  "tables"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "second_vax_period", "data", "second_vax_period_dates.rds")) 

# outcomes <- readr::read_rds(
#   here::here("output", "lib", "outcomes.rds")
# )
outcomes <- c("covidadmitted", "covidemergency")
names(outcomes) <- c("from APCS", "from ECDS")

subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds")
)

################################################################################
# function for printing glance
print_table <- function(outcome) {
  
  model_glance <- bind_rows(
    lapply(
      1:6,
      function(kk)
      readr::read_rds(
        here::here("output", "models_cox", "data", glue("modelcox_glance_{comparison}_{subgroup_label}_{outcome}_{kk}.rds"))) %>%
        mutate(model = str_c(kk, model, sep = "-"))
    )
  )
    
  
  model_glance %>%
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
  
  data <-  bind_rows(
    lapply(
      1:6,
      function(kk)
        readr::read_rds(
          here::here("output", "models_cox", "data", glue("modelcox_tidy_{comparison}_{subgroup_label}_{outcome}_{kk}.rds"))) %>%
        mutate(model = str_c(kk, model, sep = "-"))
    )
  )
    
  data %>%
    filter(is.na(reference_row) | !reference_row) %>% 
    mutate(
      lower = estimate - qnorm(0.975)*robust.se,
      upper = estimate + qnorm(0.975)*robust.se
      ) %>%
    mutate(across(c(estimate, lower, upper), 
                  ~round(exp(.x), 2))) %>%
    mutate(value = str_c(estimate, " (", lower, "-", upper, ")")) %>%
    select(term, value, model) %>% 
    # mutate(across(model, 
    #               factor,
    #               levels = as.character(1:2),
    #               labels = c("unadjusted",
    #                          "demographic + clinical"))) %>%
    pivot_wider(names_from = "model", values_from = "value") %>%
    kableExtra::kable("pipe", caption = outcome)
  
}

capture.output(
  map(outcomes, function(x) try(print_coefficients(outcome = x))),
  file = here::here("output", "models_cox", "tables", glue("modelcox_coefficients_{comparison}_{subgroup_label}.txt")),
  append=FALSE
)


