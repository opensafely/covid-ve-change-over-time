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

# study_parameters <- readr::read_rds(here::here("output", "lib", "study_parameters.rds"))

data_comparisons <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds"))

data_outcomes <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_outcomes.rds"))

# combine outcomes 
# TODO
# data_comparisons_combined <- data_comparisons %>%
#   mutate(
#     postest_date = if_else(
#       !is.na(postest_date) ~ postest_date,
#       !is.na(covidadmitted_date) ~ covidadmitted_date - days(7),
#       !is.na(coviddeath_date) ~ coviddeath_date - days(14)
#     ))

for (b in unique(data_comparisons$brand)) {
  
  cat(rep("-",40), "\n")
  cat(glue("Brand = {b}:"), "\n")
  
  # derive data_tte  
  data_tte <- data_comparisons %>%
    filter(brand %in% b) %>%
    select(patient_id, comparison, arm, time_zero_date, end_fu_date) %>%
    left_join(data_outcomes %>% 
                select(patient_id, 
                       starts_with(outcome), 
                       dereg_date, noncoviddeath_date), 
              by = "patient_id") %>%
    mutate(across(c(starts_with(outcome),dereg_date, noncoviddeath_date),
                  ~ if_else(
                    !is.na(.x) & time_zero_date < .x & .x <= end_fu_date,
                    .x,
                    as.Date(NA_character_)
                  ))) %>%
    group_by(comparison) %>%
    mutate(across(ends_with("date"),
                  ~ as.integer(.x - min(time_zero_date)))) %>%
    rename_at(vars(ends_with("_date")),
              ~ str_remove(.x, "_date")) %>%
    mutate(
      tte = pmin(!! sym(outcome), dereg, noncoviddeath, end_fu, na.rm = TRUE),
      status = if_else(
        !is.na(!! sym(outcome)) & !! sym(outcome) == tte,
        TRUE,
        FALSE
      )) %>%
    ungroup() %>%
    select(patient_id, arm, comparison, tstart = time_zero, tstop = tte, status) 
  
  counts <- data_tte %>%
    group_by(comparison) %>%
    summarise(
      n_rows = n(),
      n_events = sum(status)) %>%
    ungroup() %>%
    mutate(events_per_1000 = 100*n_events/n_rows)
  
  # checks
  cat(" \n")
  cat("summary of data_tte:", "\n")
  print(counts)
  cat("\n", glue("memory usage = ", format(object.size(data_tte), units="MB", standard="SI", digits=3L)), "\n")
  
  stopifnot("tstart should be  >= 0 in data_tte" = data_tte$tstart>=0)
  stopifnot("tstop - tstart should be strictly > 0 in data_tte" = data_tte$tstop - data_tte$tstart > 0)
  
  readr::write_rds(
    data_tte,
    here::here("output", glue("jcvi_group_{group}"), "data", glue("data_tte_{b}_{outcome}.rds")),
    compress = "gz")
  
}

