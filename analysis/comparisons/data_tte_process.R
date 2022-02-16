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

arm1 <- if_else(comparison =="ChAdOx", "ChAdOx", "BNT162b2")
arm2 <- if_else(comparison == "both", "ChAdOx", "unvax")

################################################################################
# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))

# covariates data
data_covariates <- readr::read_rds(
  here::here("output", "data", "data_covariates.rds")) %>%
  filter(arm %in% c(arm1, arm2))

# processed data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

data_all <- data_covariates %>%
  left_join(
    data_processed,
    by = "patient_id"
  )

# # outcomes for exclusions prior to comparison 1
# data_prior_outcomes <- readr::read_rds(
#   here::here("output", "data", "data_processed.rds")) %>%
#   select(patient_id, starts_with(unname(outcomes))) %>%
#   rename_with(~ glue("prior_{.x}"), starts_with(unname(outcomes)))

# # read any test data
# data_tests <- readr::read_rds(
#   here::here("output", "data", "data_tests.rds")) %>%
#   select(patient_id, matches("any_test_\\d\\_date")) %>%
#   pivot_longer(
#     cols = -patient_id,
#     names_pattern = "^(.*)_(\\d+)_date",
#     names_to = c(NA, "comparison"),
#     values_to = "anytest_date",
#     values_drop_na = TRUE
#   )

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
if ("ChAdOx" %in% c(arm1, arm2)) {
  
  
  
}

subgroups <- c(subgroups, "all")

fs::dir_create(here::here("output", "tte", "data"))
fs::dir_create(here::here("output", "tte", "tables"))

################################################################################


################################################################################

data_comparisons <- data_all %>%
  select(patient_id, k, arm, subgroup, split,
         start_k_date, end_k_date, dereg_date,
         all_of(str_c(outcomes, "_date"))) 


subgroups_1 <- data_all %>%
  filter(arm %in% arm1) %>%
  distinct(subgroup) %>%
  unlist() %>% unname() %>% as.character()
subgroups_2 <- data_all %>%
  filter(arm %in% arm2) %>%
  distinct(subgroup) %>%
  unlist() %>% unname() %>% as.character()
subgroups <- c(intersect(subgroups_1, subgroups_2))

data_comparisons <- local({
  
  data_arm1 <-  data_covariates %>%
    select(patient_id, comparison, arm, subgroup, start_k_date, end_k_date, anytest_date, split)
  
  data_arm2 <-  readr::read_rds(
    here::here("output", "comparisons", "data", glue("data_comparisons_{arm2}.rds"))) %>%
    select(patient_id, comparison, arm, subgroup, start_fu_date, end_fu_date,
           dereg_date, death_date,
           starts_with(unname(outcomes)))
  
  subgroups_1 <- unique(as.character(data_arm1$subgroup))
  subgroups_2 <- unique(as.character(data_arm2$subgroup))
  subgroups <- c(intersect(subgroups_1, subgroups_2), "all")
  
  bind_rows(data_arm1, data_arm2) %>%
    filter(subgroup %in% subgroups)

})

data <- data_comparisons %>%
  left_join(data_prior_outcomes, by = "patient_id") %>%
  left_join(data_tests, by = c("patient_id", "comparison"))

################################################################################
# generates and saves data_tte and tabulates event counts 
# returns tables of events
derive_data_tte <- function(
  .data, 
  outcome
  ) {
  
  # remove comparisons for which outcome has occurred before the patient's first comparison
  # (if outcome is anytest, only exclude if previous postest)
  outcome_exclude <- if_else(
    outcome == "anytest",
    "postest",
    outcome)
  
  data_inc <- .data %>%
    filter(
      # allow for the fact that the first comparison for the even unvax arm is 2
      (!(arm %in% "unvax") & comparison %in% "1") |
        (arm %in% "unvax" & comparison %in% c("1", "2"))  
    ) %>%
    filter(
      # remove people who have experienced the outcome before first comparison
      is.na(!! sym(glue("prior_{outcome_exclude}_date"))) |
        start_fu_date < !! sym(glue("prior_{outcome_exclude}_date"))
    ) %>%
    select(patient_id)
  
  # keep the selected patients from .data
  data_tte_0 <- data_inc %>%
    left_join(.data, by = "patient_id") %>%
    select(patient_id, comparison, arm, subgroup, start_fu_date, end_fu_date, 
           dereg_date, death_date, # for censoring
           # when outcome is anytest, keep postest too
           all_of(glue("{unique(c(outcome, outcome_exclude))}_date"))) %>%
    arrange(patient_id, comparison) 
  
  data_tte_1 <- data_tte_0 %>%
    group_by(patient_id) %>%
    # remove comparisons for which outcome_exclude has occurred before start_fu_date
    mutate(
      event_seq = cumsum(cumsum(!is.na(!! sym(glue("{outcome_exclude}_date")))))
    ) %>%
    ungroup() %>%
    filter(event_seq <= 1)
  
  data_tte_2 <- data_tte_1 %>%
    # new time-scale: time since earliest start_fu_date in data
    mutate(across(ends_with("_date"),
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
  stopifnot("tstart should be  >= 0 in data_tte_2" = data_tte_2$tstart>=0)
  stopifnot("tstop - tstart should be strictly > 0 in data_tte_2" = data_tte_2$tstop - data_tte_2$tstart > 0)
  
  # subgroups in .data
  subgroup <- unique(as.character(.data$subgroup))
  if (length(subgroup) > 1) subgroup <- "all"
  subgroup_label <- which(subgroups == subgroup)
  
  # save data_tte
  readr::write_rds(
    data_tte_2,
    here::here("output", "tte", "data", glue("data_tte_{comparison}_{subgroup_label}_{outcome}.rds")),
    compress = "gz")
  
  # tabulate events per comparison and save
  table_events <- data_tte_2 %>%
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
  


