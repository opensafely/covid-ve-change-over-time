

################################################################################

library(tidyverse)
library(survival)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobs <- FALSE
  group <- "02"
  
} else{
  removeobs <- TRUE
  group <- args[[1]]
}

study_parameters <- readr::read_rds(here::here("output", "lib", "study_parameters.rds"))

data_comparisons <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_comparisons.rds"))

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
  
  outcomes <- c("postest", "covidadmitted", "coviddeath", "death")
  
  n_outcomes <- data_comparisons %>%
    filter(brand %in% b) %>%
    summarise(across(all_of(str_c(outcomes, "_date")),
                     ~ sum(!is.na(.x)))) %>%
    pivot_longer(cols = everything()) 
  
  cat("Number of individuals with each outcome:", "\n")
  n_outcomes %>% 
    mutate(include = value > study_parameters$outcome_threshold) %>%
    # round to nearest 10
    mutate(across(value, ~round(.x,-1))) %>% 
    rename(n=value) %>%
    print()
  
  
  
}










outcomes <- c("postest", "covidadmitted", "coviddeath", "death")censor <- c("noncoviddeath", "dereg")


dates <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), 1)

test <- data_comparisons %>%
  filter(comparison==1, brand == "BNT162b2") %>%
  mutate(origin = min(time_zero_date)) %>%
  mutate(
    postest_date = sample(dates, size=nrow(.), replace=TRUE),
    covidadmitted_date = sample(dates, size=nrow(.), replace=TRUE),
    coviddeath_date = sample(dates, size=nrow(.), replace=TRUE),
    death_date = sample(dates, size=nrow(.), replace=TRUE)) %>%
  mutate(across(c(all_of(str_c(c(outcomes, censor), "_date"))),
                ~ if_else(
                  .x <= time_zero_date | .x > end_fu_date,
                  NA_integer_,
                  as.integer(.x - origin)))) %>%
  mutate(across(c(time_zero_date, end_fu_date),
                ~ as.integer(.x - origin))) 

  mutate(across(c(time_zero_date, end_fu_date, 
                  all_of(str_c(c(outcomes, censor), "_date"))),
                ~ as.integer(.x - origin)))

library(survival)
test_tmerge <- as_tibble(
  tmerge(
    
    data1 = test %>% select(patient_id),
    data2 = test,
    id = patient_id,
    
    tstart = time_zero_date,
    tstop = end_fu_date,
    
    postest = event(postest_date),
    covidadmitted = event(covidadmitted_date),
    coviddeath = event(coviddeath_date),
    death = event(death_date)
    )
  )


tte <- function(.data, var_string) {
  
  name <- str_remove(var_string, "_date")
  
  .data %>% 
    rename(temp = var_string) %>%
    mutate(
      !! glue("{name}_status") := case_when(
        is.na(temp) ~ 0L,
        time_zero_date < temp & temp <= end_fu_date ~ 1L,
        TRUE ~ 0L),
      !! glue("{name}_tte") := case_when(
        is.na(temp) ~ as.integer(end_fu_date - origin),
        time_zero_date < temp & temp <= end_fu_date ~ as.integer(temp - origin),
        TRUE ~ as.integer(end_fu_date - origin))) %>%
    select(-temp)

  }

data_tte <- data_covs %>%
  select(patient_id, 
         elig_date, region, brand, arm, time_zero_date, end_fu_date, comparison,
         coviddeath_date, noncoviddeath_date, death_date, dereg_date) %>%
  left_join(
    readr::read_rds(
      here::here("output", glue("jcvi_group_{group}"),  "data", "data_long_postest_dates.rds")
      ) %>% 
      select(patient_id, postest_date = date),
    by = "patient_id"
  ) %>%
  left_join(
    readr::read_rds(
      here::here("output", glue("jcvi_group_{group}"),  "data", "data_long_covidadmitted_dates.rds")
    ) %>% 
      select(patient_id, covidadmitted_date = date),
    by = "patient_id"
  ) %>%
  # derive origin for each brand and k
  group_by(brand, comparison) %>%
  mutate(origin = min(time_zero_date)) %>%
  ungroup() %>%
  # time to event for all outcomes and censoring events
  tte("postest_date") %>%
  tte("covidadmitted_date") %>%
  tte("coviddeath_date") %>%
  tte("death_date") %>%
  tte("noncoviddeath_date") %>%
  tte("dereg_date") %>%
  # convert time_zero and end_fu to days since origin
  mutate(
    time_zero = as.integer(time_zero_date - origin),
    end_fu = as.integer(end_fu_date - origin),
    ) 

readr::write_rds(
  data_tte, 
  here::here("output", glue("jcvi_group_{group}"),  "data", "data_tte.rds"), 
  compress="gz")
