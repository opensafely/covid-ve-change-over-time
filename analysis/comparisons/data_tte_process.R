################################################################################

# This script:
# creates time-to-event data for the given outcome

################################################################################

library(tidyverse)
library(glue)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  group <- "02"
  outcome <- "postest"
  
} else{
  group <- args[[1]]
  outcome <- args[[2]]
}

fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "tables"))

data_comparisons <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds"))

data_outcomes <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_outcomes.rds"))

median_times_between_outcomes <- readr::read_csv(
  here::here("output", glue("jcvi_group_{group}"), "data", "median_times_between_outcomes.csv"))

# combine outcomes 
# for patients that have covidadmitted but not postest (and coviddeath but not covidadmitted or postest)
# impute the median time between outcomes in those patients who do have the upstream outcomes
data_outcomes_combined <- data_outcomes %>%
  mutate(across(postest_date,
                ~ case_when(
                  !is.na(.x) ~ .x,
                  !is.na(covidadmitted_date) ~ covidadmitted_date + median_times_between_outcomes[median_times_between_outcomes$name == "postest and covidadmitted",]$median,
                  !is.na(coviddeath_date) ~ coviddeath_date + median_times_between_outcomes[median_times_between_outcomes$name == "postest and coviddeath",]$median,
                  TRUE ~ .x            
                ))) %>%
  mutate(across(covidadmitted_date,
                ~ case_when(
                  !is.na(.x) ~ .x,
                  !is.na(coviddeath_date) ~ coviddeath_date + median_times_between_outcomes[median_times_between_outcomes$name == "covidadmitted and coviddeath",]$median,
                  TRUE ~ .x            
                )))

for (b in unique(data_comparisons$brand)) {
  
  cat(rep("-",40), "\n")
  cat(glue("Brand = {b}:"), "\n")
  
  # derive data_tte  
  data_tte <- data_comparisons %>%
    filter(brand %in% b) %>%
    select(patient_id, comparison, arm, start_fu_date, end_fu_date) %>%
    left_join(data_outcomes_combined %>% 
                select(patient_id, 
                       starts_with(outcome), 
                       dereg_date, noncoviddeath_date), 
              by = "patient_id") %>%
    mutate(across(c(starts_with(outcome), dereg_date, noncoviddeath_date),
                  ~ if_else(
                    !is.na(.x) & start_fu_date < .x & .x <= end_fu_date,
                    .x,
                    as.Date(NA_character_)
                  ))) %>%
    group_by(comparison) %>%
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
    ungroup() %>%
    select(patient_id, arm, comparison, tstart = start_fu, tstop = tte, status) 
  
  # checks
  cat(" \n")
  cat("\n", glue("memory usage = ", format(object.size(data_tte), units="MB", standard="SI", digits=3L)), "\n")
  
  stopifnot("tstart should be  >= 0 in data_tte" = data_tte$tstart>=0)
  stopifnot("tstop - tstart should be strictly > 0 in data_tte" = data_tte$tstop - data_tte$tstart > 0)
  
  events_per_personyears <- function(.data) {
    .data %>%
      mutate(days = tstop-tstart) %>%
      group_by(comparison, arm) %>%
      summarise(
        n = n(),
        personyears = sum(days)/365,
        events = sum(status),
        .groups = "keep"
      ) %>%
      ungroup() %>%
      mutate(across(c(n, events, personyears), 
                    # round to nearest 10
                    ~ scales::comma(round(.x, -1), accuracy = 1))) %>%
      mutate(value = str_c(events, " / ", personyears, " (",n,")")) %>%
      select(comparison, arm, value) %>%
      pivot_wider(names_from = arm, values_from = value) %>%
      kableExtra::kable(
        "pipe",
        caption = glue("{outcome} events / person-years (n)")
      )
  }
  
  capture.output(
    data_tte %>% events_per_personyears,
    file = here::here("output", glue("jcvi_group_{group}"), "tables", glue("{b}_{outcome}_incidence.txt")),
    append=FALSE
  )
  
  readr::write_rds(
    data_tte,
    here::here("output", glue("jcvi_group_{group}"), "data", glue("data_tte_{b}_{outcome}.rds")),
    compress = "gz")
  
}

