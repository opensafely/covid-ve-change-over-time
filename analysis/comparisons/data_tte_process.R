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

################################################################################
# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroups <- c(subgroups, "all")

fs::dir_create(here::here("output", "tte", "data"))
fs::dir_create(here::here("output", "tte", "tables"))

################################################################################
arm1 <- if_else(comparison =="ChAdOx", "ChAdOx", "BNT162b2")
arm2 <- if_else(comparison == "both", "ChAdOx", "unvax")

################################################################################

derive_data <- function(
  arm_1,
  arm_2
) {
  
  data_arm1 <-  readr::read_rds(
    here::here("output", "comparisons", "data", glue("data_comparisons_{arm1}.rds"))) %>%
    select(patient_id, comparison, arm, subgroup, start_fu_date, end_fu_date,
           dereg_date, death_date,
           all_of(str_c(outcomes, "_date")))
  
  
  data_arm2 <-  readr::read_rds(
    here::here("output", "comparisons", "data", glue("data_comparisons_{arm2}.rds"))) %>%
    select(patient_id, comparison, arm, subgroup, start_fu_date, end_fu_date,
           dereg_date, death_date,
           all_of(str_c(outcomes, "_date")))
  
  subgroups_1 <- unique(as.character(data_arm1$subgroup))
  subgroups_2 <- unique(as.character(data_arm2$subgroup))
  subgroups <- c(intersect(subgroups_1, subgroups_2), "all")
  
  data <- bind_rows(data_arm1, data_arm2) %>%
    filter(subgroup %in% subgroups)

}

data <- derive_data(arm1, anm2)

################################################################################
# generates and saves data_tte and tabulates event counts 
# returns tables of events
derive_data_tte <- function(
  .data, 
  outcome
  ) {
  
  # subgroups in .data
  subgroup <- unique(as.character(.data$subgroup))
  if (length(subgroup) > 1) subgroup <- "all"
  subgroup_label <- which(subgroups == subgroup)
  
  # derive data_tte
  data_tte <- .data %>%
    select(patient_id, comparison, arm, subgroup, start_fu_date, end_fu_date, 
           dereg_date, death_date, # for censoring
           matches(str_c(outcome, "_date"))) %>%
    arrange(patient_id, comparison) %>%
    group_by(patient_id) %>%
    # remove comparisons for which outcome has occurred before start_fu_date
    mutate(
      event_seq = cumsum(cumsum(!is.na(!! sym(str_c(outcome, "_date")))))
      ) %>%
    ungroup() %>%
    filter(
      event_seq <= 1
    ) %>%
    # new time-scale: time since earliest start_fu_date in data
    mutate(across(ends_with("date"),
                  ~ as.integer(.x - min(start_fu_date)))) %>%
    rename_at(vars(ends_with("_date")),
              ~ str_remove(.x, "_date")) %>%
    mutate(
      tte = pmin(!! sym(outcome), dereg, death, end_fu, na.rm = TRUE),
      status = if_else(
        !is.na(!! sym(outcome)) & !! sym(outcome) == tte,
        TRUE,
        FALSE
      )) %>%
    select(patient_id, arm, comparison, tstart = start_fu, tstop = tte, status) %>%
    arrange(patient_id, comparison) 
  
  # checks
  stopifnot("tstart should be  >= 0 in data_tte" = data_tte$tstart>=0)
  stopifnot("tstop - tstart should be strictly > 0 in data_tte" = data_tte$tstop - data_tte$tstart > 0)
  
  # save data_tte
  readr::write_rds(
    data_tte,
    here::here("output", "tte", "data", glue("data_tte_{comparison}_{subgroup_label}_{outcome}.rds")),
    compress = "gz")
  
  # tabulate events per comparison and save
  table_events <- data_tte %>%
    mutate(days = tstop-tstart) %>%
    group_by(comparison, arm) %>%
    summarise(
      n = n(),
      personyears = sum(days)/365.25,
      events = sum(status),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    mutate(across(c(n, events, personyears),
                  # round to nearest 10
                  ~ scales::comma(round(.x, -1), accuracy = 1))) %>%
    mutate(value = str_c(events, " / ", personyears)) %>%
    select(comparison, arm, value) %>%
    rename(k = comparison) %>%
    pivot_wider(names_from = arm, values_from = value) %>%
    mutate(
      subgroup = subgroup,
      outcome = outcome
    )
  
  return(table_events)
  
}

################################################################################
# apply derive_data_tte for all comparisons, and both for all subgroups and split by subgroup

table_events <- 
  lapply(
    splice(data, as.list(data %>% group_split(subgroup))),
    function(y)
      lapply(
        outcomes,
        function(z)
          try(y %>% derive_data_tte(outcome = z))
      )
  )

readr::write_rds(
  table_events,
  here::here("output", "tte", "tables", glue("event_counts_{comparison}.rds")),
  compress = "gz")
