

################################################################################

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

vax_brand <- "BNT162b2"

outcomes <- c("postest", "covidadmitted", "coviddeath", "death")
censor <- c("noncoviddeath", "dereg")

data_covs <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_covs.rds"))

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
