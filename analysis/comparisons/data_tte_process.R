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
  comparison <- "ChAdOx"
  
} else{
  comparison <- args[[1]]
}

################################################################################
# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))

# read any test data
data_tests_1 <- readr::read_rds(
  here::here("output", "data", "data_tests.rds")) %>%
  select(patient_id, matches("any_test_\\d\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_pattern = "^(.*)_(\\d+)_date",
    names_to = c(NA, "comparison"),
    values_to = "anytest_date",
    values_drop_na = TRUE
  )

if (comparison != "both") {
  
  data_tests_unvax <- bind_rows(
    data_tests_1 %>%
    mutate(across(comparison,
                  ~ case_when(
                    .x %in% c("1","2") ~ "1",
                    .x %in% c("3","4") ~ "3",
                    .x %in% c("5", "6") ~ "5",
                    TRUE ~ NA_character_
                  ))),
    data_tests_1 %>%
      mutate(across(comparison,
                    ~ case_when(
                      .x %in% c("2","3") ~ "2",
                      .x %in% c("4","5") ~ "4",
                      .x %in% c("6", "7") ~ "6",
                      TRUE ~ NA_character_
                    )))
  ) %>%
    filter(!is.na(comparison)) %>%
    group_by(patient_id, comparison) %>%
    summarise(anytest_date = min(anytest_date), .groups = "keep") %>%
    ungroup()
  
}

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
    left_join(data_tests_1, by = c("patient_id", "comparison")) %>%
    select(patient_id, comparison, arm, subgroup, start_fu_date, end_fu_date,
           dereg_date, death_date,
           all_of(str_c(outcomes, "_date")))
  
  if (arm2 == "unvax") {
    data_tests_2 <- data_tests_unvax
  } else {
    data_tests_2 <- data_tests_1
  }
  
  data_arm2 <-  readr::read_rds(
    here::here("output", "comparisons", "data", glue("data_comparisons_{arm2}.rds"))) %>%
    left_join(data_tests_2, by = c("patient_id", "comparison")) %>%
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
    group_by(comparison, arm) %>%
    summarise(
      n = n(),
      events = sum(status),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    mutate(outcome = outcome,
           subgroup = subgroup)
  
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

table_events <- bind_rows(
  unlist(table_events, recursive = FALSE)
)

readr::write_rds(
  table_events,
  here::here("output", "tte", "tables", glue("event_counts_{comparison}.rds")),
  compress = "gz")
  


