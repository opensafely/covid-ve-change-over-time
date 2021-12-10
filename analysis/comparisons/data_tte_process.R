

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

outcomes <- c("postest", "covidadmitted", "coviddeath", "noncoviddeath", "death")

data_covariates <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_covariates.rds"))

tte <- function(.data, var_string) {
  
  name <- str_remove(var_string, "_date")
  
  .data %>% 
    rename(temp = var_string) %>%
    mutate(
      !! glue("{name}_ind") := case_when(
        is.na(temp) ~ 0L,
        time_zero < temp & temp <= end_fu_date ~ 1L,
        TRUE ~ 0L),
      !! glue("{name}_tte") := case_when(
        is.na(temp) ~ as.integer(end_fu_date - origin),
        time_zero < temp & temp <= end_fu_date ~ as.integer(temp - origin),
        TRUE ~ as.integer(end_fu_date - origin))) %>%
    select(-temp)

  }

data_tte <- data_covariates %>%
  select(patient_id, 
         elig_date, region, brand, arm, time_zero, end_fu_date, k,
         coviddeath_date, death_date) %>%
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
  group_by(brand, k) %>%
  mutate(origin = min(time_zero)) %>%
  ungroup() %>%
  tte("postest_date")


  
  
  filter(brand %in% vax_brand) %>%
  left_join(input_covs %>%
              select(patient_id, positive_test_date),
            by = "patient_id") %>%
  mutate(across(positive_test_date,
                ~ if_else(.x <= time_zero | .x > end_fu_date,
                          NA_Date_,
                          .x))) %>%
  group_by(brand, k) %>%
  mutate(origin = min(time_zero)) %>%
  ungroup() %>%
  mutate(start = as.integer(time_zero - origin),
         end = if_else(
           is.na(positive_test_date),
           as.integer(end_fu_date - origin),
           as.integer(positive_test_date - origin)),
         status = as.integer(!is.na(positive_test_date)))