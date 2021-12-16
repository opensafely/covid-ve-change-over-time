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

fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "models"))

data_comparisons <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds"))

source(here::here("analysis", "lib", "model_cox.R"))

model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

################################################################################
# apply model

# for each brand for the given jcvi_group
for (b in unique(data_comparisons$brand)) {
 
  # read data_tte
  data_tte <- readr::read_rds(
    here::here("output", glue("jcvi_group_{group}"), "data", glue("data_tte_{b}_{outcome}.rds")))
  
  # filter to given brand
  data_cox <- data_tte %>%
    left_join(data_comparisons %>%
                filter(brand %in% b) %>%
                select(patient_id, comparison, elig_date, region, 
                       unname(unlist(model_varlist))),
              by = c("patient_id", "comparison")) 
  
  # apply coxph model  
  model_output <- lapply(
    0:2,
    function(x)
      try(cox_model(
        number = x, 
        filename_prefix = glue("{b}_{outcome}"))))
  
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
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_glance.csv"))) 
  
  
  model_tidy <- bind_rows(
    lapply(
      # only bind tibbles (to avoid errors in case some models did not converge)
      seq_along(model_output)[sapply(model_output, function(x) is_tibble(x[[1]]))],
      # select tidy
      function(x) model_output[[x]]$tidy
    )) %>%
    mutate(outcome = outcome)
  
  readr::write_rds(
    model_tidy,
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_tidy.rds")), 
    compress="gz")
  readr::write_csv(
    model_tidy,
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_tidy.csv")))
  

}

