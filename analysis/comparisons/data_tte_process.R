################################################################################

# This script:
# creates time-to-event data for the given outcome

################################################################################

library(tidyverse)
library(glue)
# library(fastDummies)

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))

fs::dir_create(here::here("output", "tte", "data"))
fs::dir_create(here::here("output", "tte", "tables"))

################################################################################
load_data <- function(arm) {
  
  readr::read_rds(
    here::here("output", "data", glue("data_comparisons_{arm}.rds"))) %>%
    select(patient_id, comparison, arm, subgroup, start_fu_date, end_fu_date) %>%
    left_join(
      readr::read_rds(
        here::here("output", "data", glue("data_outcomes_{arm}.rds"))) %>% 
        select(patient_id, all_of(str_c(outcomes, "_date")), 
               dereg_date), 
      by = "patient_id")
  
}

################################################################################
# load comparisons and outcomes data
data_BNT162b2 <- load_data("BNT162b2")

data_ChAdOx <- load_data("ChAdOx")

data_unvax <- load_data("unvax")

################################################################################

derive_data <- function(
  data_arm_1,
  data_arm_2
) {
  
  subgroups_1 <- unique(as.character(data_arm_1$subgroup))
  subgroups_2 <- unique(as.character(data_arm_2$subgroup))
  subgroups <- c(intersect(subgroups_1, subgroups_2), "all")
  
  data <- bind_rows(
    data_arm_1,
    data_arm_2
  )  %>%
    filter(subgroup %in% subgroups)

}

data_BNT162b2_unvax <- derive_data(data_BNT162b2, data_unvax)
data_ChAdOx_unvax <- derive_data(data_ChAdOx, data_unvax)
data_BNT162b2_ChAdOx <- derive_data(data_BNT162b2, data_ChAdOx)

rm(data_BNT162b2, data_ChAdOx, data_unvax)

################################################################################
# generates and saves data_tte and tabulates event counts (does not return anything)
derive_data_tte <- function(
  .data, 
  outcome
  ) {
  
  # comparison arms
  arms <- unique(as.character(.data$arm))
  if ("unvax" %in% arms) {
    comparison <- arms[which(arms != "unvax")]
  } else {
    comparison <- "both"
  }
  
  # subgroups in .data
  subgroup <- unique(as.character(.data$subgroup))
  if (length(subgroup) > 1) subgroup <- "all"
  
  # derive data_tte
  data_tte <- .data %>%
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
  stopifnot("tstart should be  >= 0 in data_tte" = data_tte$tstart>=0)
  stopifnot("tstop - tstart should be strictly > 0 in data_tte" = data_tte$tstop - data_tte$tstart > 0)
  
  # save data_tte
  readr::write_rds(
    data_tte,
    here::here("output", "tte", "data", glue("data_tte_{comparison}_{subgroup}_{outcome}.rds")),
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
    pivot_wider(names_from = arm, values_from = value) 
  
  readr::write_rds(
    table_events,
    here::here("output", "tte", "tables", glue("table_events_{comparison}_{subgroup}_{outcome}.rds")),
    compress = "gz")
  
  return(cat(glue("data_tte for comparison {comparison}, subgroup {subgroup} and outcome {outcome} completed"), "\n"))
  
}

################################################################################
# apply derive_data_tte for all comparisons, and both for all subgroups and split by subgroup

stopreturn <- lapply(
  list(
    data_BNT162b2_unvax,
    data_ChAdOx_unvax,
    data_BNT162b2_ChAdOx
  ),
  function(x) 
    lapply(
      splice(x, as.list(x %>% group_split(subgroup))),
      function(y)
        lapply(
          outcomes,
          function(z)
            try(y %>% derive_data_tte(outcome = z))
        )
    )
)



