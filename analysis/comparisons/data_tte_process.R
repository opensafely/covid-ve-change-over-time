################################################################################

# This script:
# creates time-to-event data for the given outcome

################################################################################

library(tidyverse)
library(glue)
library(fastDummies)

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

if (comparison %in% c("BNT162b2", "ChAdOx")) {
  
  trt <- str_c(comparison, "_vax")
  
  data <- bind_rows(
    
    lapply(
      c("vax", "unvax"),
      function(x)
        readr::read_rds(
          here::here("output", "data", glue("data_comparisons_{comparison}_{x}.rds"))) %>%
        select(patient_id, comparison, arm, start_fu_date, end_fu_date) %>%
        left_join(
          readr::read_rds(
            here::here("output", "data", glue("data_outcomes_{comparison}_{x}.rds"))) %>% 
            select(patient_id, 
                   starts_with(outcome), 
                   dereg_date, noncoviddeath_date), 
          by = "patient_id")
    )
  )
    
} else if (comparison %in% "both") {
  
  trt <- "BNT162b2_vax"
  
  data <- bind_rows(
    
    lapply(
      c("BNT162b2", "ChAdOx"),
      function(x)
        readr::read_rds(
          here::here("output", "data", glue("data_comparisons_{x}_vax.rds"))) %>%
        # remove certain groups for brands comparison
        filter(!(jcvi_group %in% c("01", "11", "12"))) %>%
        select(patient_id, comparison, arm, start_fu_date, end_fu_date) %>%
        left_join(
          readr::read_rds(
            here::here("output", "data", glue("data_outcomes_{x}_vax.rds"))) %>% 
            select(patient_id, 
                   starts_with(outcome), 
                   dereg_date, noncoviddeath_date), 
          by = "patient_id")
    )
  )
  
}

################################################################################
# derive data_tte ----
data_tte <- data %>%
  mutate(across(c(starts_with(outcome), dereg_date, noncoviddeath_date),
                ~ if_else(
                  !is.na(.x) & (start_fu_date < .x) & (.x <= end_fu_date),
                  .x,
                  as.Date(NA_character_)
                ))) %>%
  mutate(across(ends_with("date"),
                ~ as.integer(.x - min(start_fu_date)))) %>%
  rename_at(vars(ends_with("_date")),
            ~ str_remove(.x, "_date")) %>%
  mutate(
    tte = pmin(!! sym(outcome), dereg, noncoviddeath, end_fu, na.rm = TRUE),
    status = if_else(
      !is.na(!! sym(outcome)) & !! sym(outcome) == tte,
      TRUE,
      FALSE
    )) %>%
  select(patient_id, arm, comparison, tstart = start_fu, tstop = tte, status) %>%
  arrange(patient_id, comparison) 

# checks
cat(" \n")
cat("\n", glue("memory usage = ", format(object.size(data_tte), units="MB", standard="SI", digits=3L)), "\n")

stopifnot("tstart should be  >= 0 in data_tte" = data_tte$tstart>=0)
stopifnot("tstop - tstart should be strictly > 0 in data_tte" = data_tte$tstop - data_tte$tstart > 0)

readr::write_rds(
  data_tte,
  here::here("output", "data", glue("data_tte_{comparison}_{outcome}.rds")),
  compress = "gz")
