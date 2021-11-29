
# # # # # # # # # # # # # # # # # # # # #
# This script:
# creates metadata for aspects of the study design
# # # # # # # # # # # # # # # # # # # # #

# Import libraries ----
library(tidyverse)
library(lubridate)
library(glue)

# create output directories ----
fs::dir_create(here::here("output", "lib"))

# create study_parameters ----

study_parameters <-
  list(
    n_comparisons = 3, # the number of comparisons for each sequence
    ref_age_1 = "2021-03-31", # reference date for calculating age for phase 1 groups
    ref_age_2 = "2021-07-01", # reference date for calculating age for phase 2 groups
    ref_cev = "2021-01-18", # reference date for calculating eligibility for phase 1 group 4 (CEV)
    ref_ar = "2021-02-15", # reference date for calculating eligibility for phase 1 group 5 (at-risk)
    pandemic_start = "2020-01-01", # rough start date for pandemic in UK
    start_date = "2020-12-08", # start of phase 1 vaccinations
    start_date_pfizer = "2020-12-08",
    start_date_az = "2021-01-04",
    start_date_moderna = "2021-03-04",
    end_date = "2021-09-15" # last date of available vaccination data. NEED TO ALSO CHECK END DATES FOR OTHER DATA SOURCES
  )

readr::write_rds(study_parameters, here::here("output", "lib", "study_parameters.rds"))
jsonlite::write_json(study_parameters, path = here::here("output", "lib", "study_parameters.json"), auto_unbox = TRUE, pretty=TRUE)

# create jcvi_groups ----
jcvi_groups <- 
tribble(
    ~group, ~definition,
    "01", "longres_group AND age_1 > 65",
    "02", "age_1 >=80",
    "03", "age_1 >=75",
    "04", "age_1 >=70 OR (cev_group AND age_1 >=16 AND NOT preg_group)",
    "05", "age_1 >=65",
    "06", "atrisk_group AND age_1 >=16",
    "07", "age_1 >=60",
    "08", "age_1 >=55",
    "09", "age_1 >=50",
    "10", "age_2 >=40",
    "11", "age_2 >=30",
    "12", "age_2 >=18",
    "99", "DEFAULT",
)

readr::write_csv(jcvi_groups, here::here("output", "lib", "jcvi_groups.csv"))

# create elig_dates ----
elig_dates <-
tribble(
    ~date, ~description, ~jcvi_groups,
    "2020-12-08", "jcvi_group='01' OR jcvi_group='02' OR jcvi_group='03'", "01, 02, 03",
    "2021-01-18", "jcvi_group='04'", "04",
    ###
    "2021-02-15", "jcvi_group='05' OR jcvi_group='06'", "05, 06",
    ###
    "2021-02-22", "age_1 >= 64 AND age_1 < 65", "07", 
    "2021-03-01", "age_1 >= 60 AND age_1 < 64", "07",
    ###
    "2021-03-08", "age_1 >= 56 AND age_1 < 60", "08",
    "2021-03-09", "age_1 >= 55 AND age_1 < 56", "08",
    ###
    "2021-03-19", "age_1 >= 50 AND age_1 < 55", "09",
    ###
    "2021-04-13", "age_2 >= 45 AND age_1 < 50", "10",
    "2021-04-26", "age_2 >= 44 AND age_1 < 45", "10",
    "2021-04-27", "age_2 >= 42 AND age_1 < 44", "10",
    "2021-04-30", "age_2 >= 40 AND age_1 < 42", "10",
    ###
    "2021-05-13", "age_2 >= 38 AND age_2 < 40", "11",
    "2021-05-19", "age_2 >= 36 AND age_2 < 38", "11",
    "2021-05-21", "age_2 >= 34 AND age_2 < 36", "11",
    "2021-05-25", "age_2 >= 32 AND age_2 < 34", "11",
    "2021-05-26", "age_2 >= 30 AND age_2 < 32", "11",
    ###
    "2021-06-08", "age_2 >= 25 AND age_2 < 30", "12",
    "2021-06-15", "age_2 >= 23 AND age_2 < 25", "12",
    "2021-06-16", "age_2 >= 21 AND age_2 < 23", "12",
    "2021-06-18", "age_2 >= 18 AND age_2 < 21", "12",
    "2100-12-31", "DEFAULT", "NA",
)

readr::write_csv(elig_dates, here::here("output", "lib", "elig_dates.csv"))

# create start and end dates for comparison_1 ----
# this is just temporary, as these dates will be saved in the file
# output/lib/second_vax_period_dates.csv
# which is output from
# analysis/eda_index_dates/plot_2nd_vax_dates.R
comparison_1_dates <-
  tribble(
    ~elig_date, 
    "2020-12-08", 
    "2021-01-18", 
    ###
    "2021-02-15", 
    ###
    "2021-02-22", 
    "2021-03-01", 
    ###
    "2021-03-08", 
    "2021-03-09", 
    ###
    "2021-03-19", 
    ###
    "2021-04-13", 
    "2021-04-26", 
    "2021-04-27", 
    "2021-04-30", 
    ###
    "2021-05-13", 
    "2021-05-19", 
    "2021-05-21", 
    "2021-05-25", 
    "2021-05-26", 
    ###
    "2021-06-08", 
    "2021-06-15", 
    "2021-06-16", 
    "2021-06-18", 
    "2100-12-31"
  ) %>%
  mutate(
    start_1_date = as.Date(elig_date, format = "%Y-%m-%d") + weeks(6),
    end_1_date = start_1_date + weeks(6),
  ) %>%
  mutate(across(elig_date,
                ~if_else(.x == "2100-12-31",
                         "DEFAULT",
                         as.character(glue("elig_date = {elig_date}"))
                ))) %>%
  distinct(elig_date, .keep_all = TRUE) 

update_comparison_dates <- function(.data, n) {
  .data %>%
    mutate("start_{n}_date" := start_1_date + days((n-1)*28),
           "end_{n}_date" := end_1_date + days((n-1)*28))
}

comparison_dates <- comparison_1_dates
for (k in 2:study_parameters$n_comparisons) {
  comparison_dates <- comparison_dates %>% update_comparison_dates(n=k)
}

readr::write_csv(comparison_dates, here::here("output", "lib", "comparison_dates.csv"))
