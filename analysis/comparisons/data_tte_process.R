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
  comparison <- "BNT162b2"
  
} else{
  comparison <- args[[1]]
}

arm1 <- if_else(comparison =="ChAdOx1", "ChAdOx1", "BNT162b2")
arm2 <- if_else(comparison == "both", "ChAdOx1", "unvax")

################################################################################
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))
outcomes_death <- outcomes[str_detect(outcomes, "death")]

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)
if ("ChAdOx1" %in% c(arm1, arm2)) {
  select_subgroups <- subgroups[subgroups != "18-39 years"]
} else {
  select_subgroups <- subgroups
}

################################################################################
# covariates data
data_covariates <- readr::read_rds(
  here::here("output", "data", "data_covariates.rds")) %>%
  filter(arm %in% c(arm1, arm2))

# processed data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
# output directories
fs::dir_create(here::here("output", "tte", "data"))
fs::dir_create(here::here("output", "tte", "tables"))

################################################################################

data <- data_covariates %>%
  left_join(
    data_processed,
    by = "patient_id"
  ) %>%
  # filter subgroups
  filter(subgroup %in% select_subgroups) %>%
  # keep only odd unvax for odd k, equiv. for even
  filter(
    is.na(split) |
      ((k %% 2) == 0 & split == "even") |
      ((k %% 2) != 0 & split == "odd")
  ) %>%
  select(patient_id, k, arm, subgroup, split,
         start_k_date, end_k_date, subsequent_vax_date, dereg_date,
         all_of(str_c(outcomes, "_date"))) 

################################################################################
# generates and saves data_tte and tabulates event counts 
# returns tables of events
derive_data_tte <- function(
  .data, 
  outcome
  ) {
  
  # remove comparisons for which outcome has occurred before the patient's first comparison
  # (if outcome is anytest, only exclude if previous postest)
  if (outcome == "anytest") {
    outcome_exclude <- "postest"
  } else if (outcome == "covidemergency") {
    outcome_exclude <- "covidadmitted" # to ensure the same sample for the hospitalisations comparison
  } else {
    outcome_exclude <- outcome
  }
  
  # function to be applied in dplyr::filter
  occurs_after_start_date <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
  data_tte <- .data %>%
    # exclude if subsequent_vax, death, dereg or outcome_exclude occurred before start of period
    filter_at(
      vars(str_c(unique(c("subsequent_vax", "dereg", "coviddeath", "noncoviddeath", outcome_exclude)), "_date")),
      all_vars(occurs_after_start_date(cov_date = ., index_date = start_k_date))
    ) %>%
    # only keep periods for which start_k_date < end_date
    filter(
      start_k_date < as.Date(study_parameters$end_date) 
    ) %>%
    # if end_k_date > end_date, replace with end_date
    mutate(across(end_k_date,
                  ~ if_else(as.Date(study_parameters$end_date) < .x,
                            as.Date(study_parameters$end_date),
                            .x))) %>%
    # only keep dates for censoring and outcome variables between start_k_date and end_k_date
    mutate(across(all_of(str_c(unique(c("dereg", outcomes)), "_date")),
                  ~ if_else(
                    !is.na(.x) & (start_k_date < .x) & (.x <= end_k_date),
                    .x,
                    as.Date(NA_character_)
                  ))) %>%
    # new time-scale: time since earliest start_fu_date in data
    mutate(across(ends_with("_date"),
                  ~ as.integer(.x - min(start_k_date)))) %>%
    rename_at(vars(ends_with("_date")),
              ~ str_remove(.x, "_date")) %>%
    mutate(
      # censor follow-up time at first of the following:
      tte = pmin(!! sym(outcome), dereg, coviddeath, noncoviddeath, end_k, na.rm = TRUE),
      status = if_else(
        !is.na(!! sym(outcome)) & !! sym(outcome) == tte,
        TRUE,
        FALSE
      )) %>%
    select(patient_id, arm, k, tstart = start_k, tstop = tte, status) %>%
    arrange(patient_id, k) 
  
  # checks
  stopifnot("tstart should be  >= 0 in data_tte" = data_tte$tstart>=0)
  stopifnot("tstop - tstart should be strictly > 0 in data_tte" = data_tte$tstop - data_tte$tstart > 0)
  
  # subgroups in .data
  subgroup_current <- unique(as.character(.data$subgroup))
  subgroup_current_label <- subgroup_labels[subgroups == subgroup_current]
  
  # save data_tte
  readr::write_rds(
    data_tte,
    here::here("output", "tte", "data", glue("data_tte_{comparison}_{subgroup_current_label}_{outcome}.rds")),
    compress = "gz")
  
  # tabulate events per comparison and save
  table_events <- data_tte %>%
    mutate(person_days = tstop-tstart) %>%
    group_by(k, arm) %>%
    summarise(
      n = n(),
      person_years = sum(person_days)/365.25,
      events = sum(status),
      .groups = "keep"
    ) %>%
    ### REDACT SMALL NUMBERS ###
    # remove counts <=5 from n and events
    mutate(across(c(n, events), redactor2)) %>%
    # also redact 0 counts (not done in redactor2)
    mutate(across(c(n, events), 
                  ~if_else(.x==0,
                           NA_integer_, 
                           .x))) %>%
    # also redact person_years if n is redacted
    # if not redacting, round to 2 d.p.
    mutate(across(person_years, 
                  ~if_else(is.na(n),
                           NA_real_, 
                           round(.x, 2)))) %>%
    ungroup() %>%
    mutate(outcome = outcome,
           subgroup = subgroup_current_label) %>%
    select(subgroup, arm, outcome, k, n, person_years, events)
  
  return(table_events)
  
}

################################################################################
# apply derive_data_tte for all comparisons, and both for all subgroups and split by subgroup

table_events <- 
  lapply(
    as.list(data %>% group_split(subgroup)),
    function(y)
      lapply(
        outcomes,
        function(z)
          try(y %>% derive_data_tte(outcome = z))
      )
  )

table_events <- bind_rows(
  unlist(table_events, recursive = FALSE)
) 

readr::write_csv(
  table_events,
  here::here("output", "tte", "tables", glue("event_counts_{comparison}.csv")))
  


