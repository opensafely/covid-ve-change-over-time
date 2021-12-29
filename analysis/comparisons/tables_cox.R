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

fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "tables"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds")) %>%
  filter(jcvi_group %in% group, include) %>%
  distinct(brand, n_comparisons)


outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)

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
