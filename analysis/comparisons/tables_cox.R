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
  removeobs <- TRUE
  group <- args[[1]]
}

################################################################################

fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "tables"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds")) %>%
  filter(jcvi_group %in% group, include) %>%
  distinct(brand, n_comparisons)


for (b in as.character(unique(second_vax_period_dates$brand))) {
  
  model_tidy <- lapply(
    
  )
  
  model_tidy <- readr::read_rds(
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_tidy.rds"))) 
  
  model_tidy <- readr::read_rds(
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_tidy.rds"))) 
  
}