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
# read data
data_all <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) %>%
  filter(arm %in% c(arm1, arm2))

data_episodes <- readr::read_rds(
  here::here("output", "data", "data_episodes.rds"))

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))
source(here::here("analysis", "functions", "data_process_functions.R"))

################################################################################
# output directories
fs::dir_create(here::here("output", "tte", "data"))
fs::dir_create(here::here("output", "tte", "tables"))

################################################################################

data_comparisons <- data_all %>%
  # filter subgroups
  filter(subgroup %in% select_subgroups) %>%
  pivot_longer(
    cols = matches("\\w+_\\d_date"),
    names_to = c(".value", "k"),
    names_pattern = "(.*)_(.)_date"
  ) %>%
  rename_with(~str_c(.x, "_date"), .cols = all_of(c("start", "end", "anytest"))) %>%
  mutate(across(k, as.integer)) %>%
  # keep only odd unvax for odd k, equiv. for even
  filter(
    is.na(split) |
      ((k %% 2) == 0 & split == "even") |
      ((k %% 2) != 0 & split == "odd")
  ) %>%
  select(patient_id, k, arm, subgroup, sex, ends_with("date")) 

# derive time since start of most recent covid episode
data_prior_covid <- data_comparisons %>%
  left_join(
    data_episodes %>% select(patient_id, episode_start_date, episode_end_date), 
    by = "patient_id") %>%
  # only keep individuals with prior covid, and only keep comparison periods that start on or after episode_end_date
  # we are allowing the episode to have ended on start date (we don't count events that occur on start date)
  # if we didn't, then the event occuring on end date (which could be a covid hospitalisation) wouldn't
  # be counted as prior covid or an outcome.
  filter(!is.na(episode_end_date) & episode_end_date <= start_date) %>%
  # arrange so that most recent episodes first
  arrange(patient_id, desc(episode_start_date)) %>%
  # keep only the most recent for each individual
  distinct(.keep_all = TRUE) %>%
  # time in days between start of most recent covid episode and start date + 1 
  # (start date + 1 because this is the first date on which they can experience an outcome)
  mutate(prior_covid = as.numeric(start_date - episode_start_date) + 1) %>%
  select(patient_id, k, prior_covid)
  
# clean data to only keep:
# 1. all comparison periods for patients with no episodes
# 2. episodes that start before and end after comparison period start_date
# 3. episodes that start during a comparison period
data_in <- data_comparisons %>%
  left_join(
    data_episodes, 
    by = "patient_id") %>%
  filter(
      # 1.
      (is.na(episode_start_date)) |
      # 2.
      (episode_start_date <= start_date & start_date <= episode_end_date) |
      # 3.
      (start_date < episode_start_date & episode_start_date <= episode_end_date)
  ) %>%
  select(-episode_start_date, -episode_end_date) %>%
  left_join(data_prior_covid, 
            by = c("patient_id", "k")) %>%
  mutate(across(prior_covid,
                ~ fct_case_when(
                  is.na(.x) ~ "no prior covid", 
                  # .x <= 90 ~ "1-90 days", # there will be very few events here, as events within 90 days would be grouped with the previous episode
                  .x <= 180 ~ "1-180 days",
                  TRUE ~ "181+ days"
                ))) 

# 1. update prior covid at the start of each comparison period
# 2. define time since prior covid at start of svp + 2 weeks 
# go with option 1, as option 2 doesn't account for the fact that
# a patient could have an event in the first comparison period and the last comparison period 
# (this will become more of a concern if we extend follow-up to >6 comparison periods)

# NOTE: maybe I should edit data_episodes in data_covariates_process.R so that 
# postest_date for every covid episode is  episode_start_date??

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
      all_vars(occurs_after_start_date(cov_date = ., index_date = start_date))
    ) %>%
    # only keep periods for which start_date < end_date 
    filter(
      start_date < as.Date(study_parameters$end_date) 
    ) %>%
    # if end_date > study_parameters$end_date, replace with study_parameters$end_date
    mutate(across(end_date,
                  ~ if_else(as.Date(study_parameters$end_date) < .x,
                            as.Date(study_parameters$end_date),
                            .x))) %>%
    # only keep dates for censoring and outcome variables between start_date and end_date
    mutate(across(all_of(str_c(unique(c("dereg", outcomes)), "_date")),
                  ~ if_else(
                    !is.na(.x) & (start_date < .x) & (.x <= end_date),
                    .x,
                    as.Date(NA_character_)
                  ))) %>%
    # new time-scale: time since earliest start_fu_date in data
    mutate(across(ends_with("_date"),
                  ~ as.integer(.x - min(start_date)))) %>%
    rename_at(vars(ends_with("_date")),
              ~ str_remove(.x, "_date")) %>%
    mutate(
      # censor follow-up time at first of the following:
      tte = pmin(!! sym(outcome), dereg, coviddeath, noncoviddeath, end, na.rm = TRUE),
      status = if_else(
        !is.na(!! sym(outcome)) & !! sym(outcome) == tte,
        TRUE,
        FALSE
      )) %>%
    select(patient_id, arm, k, tstart = start, tstop = tte, status) %>%
    arrange(patient_id, k) 
  
  # checks
  stopifnot("tstart should be  >= 0 in data_tte" = data_tte$tstart>=0)
  stopifnot("tstop - tstart should be strictly > 0 in data_tte" = data_tte$tstop - data_tte$tstart > 0)
  
  # subgroups in .data
  subgroup_current <- unique(as.character(.data$subgroup))
  subgroup_current_label <- subgroup_labels[subgroups == subgroup_current]
  # sex in .data
  sex_current_label <- unique(as.character(.data$sex))
  if (length(sex_current_label)==1) subgroup_current_label <- glue("{subgroup_current_label}_{sex_current_label}")
  
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
      person_days = sum(person_days),
      events = sum(status),
      .groups = "keep"
    ) %>%
    # round n and events up to nearest 7 for disclosure control
    mutate(across(c(n, events), ~ceiling_any(.x, to=7))) %>%
    mutate(person_years = round(person_days/365, 0)) %>%
    ungroup() %>%
    mutate(outcome = outcome,
           subgroup = as.character(subgroup_current_label)) %>%
    select(subgroup, arm, outcome, k, n, person_years, events) 
  
  return(table_events)
  
}

################################################################################
# apply derive_data_tte for all comparisons, and both for all subgroups and split by subgroup

table_events_list <- 
  lapply(
    splice(
      as.list(data %>% group_split(subgroup)), 
      as.list(data %>% group_split(subgroup, sex))
      ),
    function(y)
      lapply(
        outcomes,
        function(z)
          try(y %>% derive_data_tte(outcome = z))
      )
  )

table_events <- bind_rows(
  unlist(table_events_list, recursive = FALSE)
) %>% 
  arrange(subgroup, outcome, k, arm) 

# save for releasing
readr::write_csv(
  table_events,
  here::here("output", "tte", "data", glue("event_counts_{comparison}.csv")))

# save for checking
capture.output(
  table_events %>%
    kableExtra::kable("pipe"),
  file = here::here("output", "tte", "tables", glue("event_counts_{comparison}.txt")),
  append = FALSE
)
