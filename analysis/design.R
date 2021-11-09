
# # # # # # # # # # # # # # # # # # # # #
# This script:
# creates metadata for aspects of the study design
# # # # # # # # # # # # # # # # # # # # #

# Import libraries ----
library('tidyverse')
library('here')

# create output directories ----
fs::dir_create(here("analysis", "lib"))

# create key dates ----

dates <-
  list(
    ref_age_1 = "2021-03-31", # reference date for calculating age for phase 1 groups
    ref_age_2 = "2021-07-01", # reference date for calculating age for phase 2 groups
    ref_cev = "2021-01-18", # reference date for calculating eligibility for phase 1 group 4 (CEV)
    ref_ar = "2021-02-15", # reference date for calculating eligibility for phase 1 group 5 (at-risk)
    start_date = "2020-12-08", # start of phase 1 vaccinations
    start_date_pfizer = "2020-12-08",
    start_date_az = "2021-01-04",
    start_date_moderna = "2021-03-04",
    end_date = "2021-09-15" # last date of available vaccination data
  )

readr::write_rds(dates, here::here("analysis", "lib", "dates.rds"))
jsonlite::write_json(dates, path = here::here("analysis", "lib", "dates.json"), auto_unbox = TRUE, pretty=TRUE)

# variable labels ----


## variable labels
# variable_labels <-
#   list(
#     vax1_type ~ "Vaccine type",
#     vax1_type_descr ~ "Vaccine type",
#     age ~ "Age",
#     ageband ~ "Age",
#     sex ~ "Sex",
#     ethnicity_combined ~ "Ethnicity",
#     imd_Q5 ~ "IMD",
#     region ~ "Region",
#     stp ~ "STP",
#     vax1_day ~ "Day of vaccination",
#     jcvi_group ~ "JCVI priority group"
#   ) %>%
#   set_names(., map_chr(., all.vars))
# 
# write_rds(variable_labels, here("analysis", "lib", "variable_labels.rds"))