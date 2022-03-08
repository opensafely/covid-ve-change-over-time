################################################################################
# This script:
# creates metadata for aspects of the study design
################################################################################

# Import libraries ----
library(tidyverse)
library(lubridate)
library(glue)

print(Sys.getenv("OPENSAFELY_BACKEND"))

# create output directories ----
fs::dir_create(here::here("output", "lib"))

################################################################################
# create study_parameters ----
study_parameters <-
  list(
    seed = 123456L,
    n = 100000L, # number of individuals in dummy data
    max_comparisons = 6L, # the number of comparisons for each sequence
    n_threshold = integer(), # the number of individuals with a second dose in the second vaccination period for a given jcvi_group and brand to include comparison
    ref_age_1 = "2021-03-31", # reference date for calculating age for phase 1 groups
    ref_age_2 = "2021-07-01", # reference date for calculating age for phase 2 groups
    ref_cev = "2021-01-18", # reference date for calculating eligibility for phase 1 group 4 (CEV)
    ref_ar = "2021-02-15", # reference date for calculating eligibility for phase 1 group 5 (at-risk)
    pandemic_start = "2020-01-01", # rough start date for pandemic in UK
    start_date = "2020-12-08", # start of phase 1 vaccinations
    start_date_pfizer = "2020-12-08",
    start_date_az = "2021-01-04",
    start_date_moderna = "2021-03-04",
    end_date = "2021-12-15" # last date of available data
  ) 

# use lower thresholds if not running in the server
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  study_parameters$n_threshold <- 100L
  # study_parameters$outcome_threshold <- 10L
} else {
  study_parameters$n_threshold <- 1000L
  # study_parameters$outcome_threshold <- 100L
}

readr::write_rds(study_parameters, here::here("output", "lib", "study_parameters.rds"))
jsonlite::write_json(study_parameters, path = here::here("output", "lib", "study_parameters.json"), auto_unbox = TRUE, pretty=TRUE)

################################################################################
# create jcvi_groups ----
jcvi_groups <- 
tribble(
    ~group, ~definition,
    "01", "longres_group AND age_1 > 65",
    "02", "age_1 >=80",
    "03", "age_1 >=75",
    "04", "age_1 >=70 OR (cev_group AND age_1 >=16)",
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

################################################################################
# create elig_dates ----
# group elig_date if within 7 days of previous elig_date (within jcvi_group)
elig_dates <-
tribble(
    ~date, ~description, ~jcvi_groups,
    "2020-12-08", "jcvi_group='01' OR jcvi_group='02' OR jcvi_group='03'", "01, 02, 03", #TODO
    "2021-01-18", "jcvi_group='04'", "04",
    ###
    "2021-02-15", "jcvi_group='05' OR jcvi_group='06'", "05, 06",
    ###
    "2021-02-22", "age_1 >= 64 AND age_1 < 65", "07", 
    "2021-03-01", "age_1 >= 60 AND age_1 < 64", "07",
    ###
    # combine 2 rows as < 7 days between elig_dates
    "2021-03-08", "age_1 >= 55 AND age_1 < 60", "08",
    # "2021-03-08", "age_1 >= 56 AND age_1 < 60", "08",
    # "2021-03-09", "age_1 >= 55 AND age_1 < 56", "08",
    ###
    "2021-03-19", "age_1 >= 50 AND age_1 < 55", "09",
    ###
    "2021-04-13", "age_2 >= 45 AND age_1 < 50", "10",
    # combine 3 rows as < 7 days between elig_dates
    "2021-04-26", "age_2 >= 40 AND age_1 < 45", "10",
    # "2021-04-26", "age_2 >= 44 AND age_1 < 45", "10",
    # "2021-04-27", "age_2 >= 42 AND age_1 < 44", "10",
    # "2021-04-30", "age_2 >= 40 AND age_1 < 42", "10",
    ###
    # combine 2 rows as < 7 days between elig_dates
    "2021-05-13", "age_2 >= 36 AND age_2 < 40", "11",
    # "2021-05-13", "age_2 >= 38 AND age_2 < 40", "11",
    # "2021-05-19", "age_2 >= 36 AND age_2 < 38", "11",
    # combine 3 rows as < 7 days between elig_dates
    "2021-05-21", "age_2 >= 30 AND age_2 < 36", "11",
    # "2021-05-21", "age_2 >= 34 AND age_2 < 36", "11",
    # "2021-05-25", "age_2 >= 32 AND age_2 < 34", "11",
    # "2021-05-26", "age_2 >= 30 AND age_2 < 32", "11",
    ###
    "2021-06-08", "age_2 >= 25 AND age_2 < 30", "12",
    # combine 3 rows as < 7 days between elig_dates
    "2021-06-15", "age_2 >= 18 AND age_2 < 25", "12",
    # "2021-06-15", "age_2 >= 23 AND age_2 < 25", "12",
    # "2021-06-16", "age_2 >= 21 AND age_2 < 23", "12",
    # "2021-06-18", "age_2 >= 18 AND age_2 < 21", "12",
    "2100-12-31", "DEFAULT", "NA",
) 

readr::write_csv(elig_dates, here::here("output", "lib", "elig_dates.csv"))

################################################################################
# create regions ----
regions <- tribble(
  ~region, ~ratio,
  "North East", 0.1,
  "North West", 0.1,
  "Yorkshire and The Humber", 0.1,
  "East Midlands", 0.1,
  "West Midlands", 0.1,
  "East", 0.1,
  "London", 0.2,
  "South West", 0.1,
  "South East", 0.1
)

readr::write_csv(regions, here::here("output", "lib", "regions.csv"))

################################################################################
# varlists for cox models

clinical <-c(
  "BMI" = "bmi",
  "Chronic respiratory disease" = "chronic_respiratory_disease", 
  "Chronic heart disease" = "chronic_heart_disease", 
  "Chronic liver disease" = "chronic_liver_disease", 
  "Chronic kidney disease" = "ckd", 
  "Chronic neurological disease" = "chronic_neuro_inc_ld",
  "Diabetes" = "diabetes",
  "Immunosuppression" = "any_immunosuppression",
  "Learning disability" = "ld_inc_ds_and_cp",
  "Serious mental illness" = "sev_ment",
  "Shielding criteria met" = "cev", 
  "Morbidity count" = "multimorb",
  "Flu vaccine in previous 5 years" = "flu_vaccine",
  "Resident in long-term residential home" = "longres", 
  "Housebound" = "housebound",
  "Number of SARS-CoV-2 tests between 2020-05-18 and min_elig_date" = "test_hist_n",
  "Pregnancy" = "pregnancy"
)

# clinical <- c(
#   "BMI" = "bmi",
#   "Heart failure" = "heart_failure", 
#   "Other heart disease" = "other_heart_disease", 
#   "Dialysis" = "dialysis",
#   "Diabetes" = "diabetes",
#   "Chronic liver disease" = "chronic_liver_disease", 
#   "COPD" = "current_copd",
#   "Other respiratory disease" = "other_respiratory", 
#   "Lung cancer" = "lung_cancer",
#   "Haematological cancer" = "haematological_cancer",
#   "Cancer excl. lung, haemo" = "cancer_excl_lung_and_haem",
#   "Immunosuppressed" = "any_immunosuppression",
#   "Dementia"  = "dementia", 
#   "Other neurological conditions" = "other_neuro_conditions",
#   "Learning disabilities" = "ld_inc_ds_and_cp",
#   "Serious mental illness" = "psychosis_schiz_bipolar",
#   "Morbidity count" = "multimorb",
#   "Shielding criteria met" = "shielded", 
#   "Flu vaccine in previous 5 years" = "flu_vaccine",
#   "Resident in long-term residential home" = "longres", 
#   "Number of SARS-CoV-2 tests between 2020-05-18 and min_elig_date" = "test_hist_1_n",
#   "Pregnancy" = "pregnancy"
# )

demographic <- c(
  "Age" = "age",
  "Sex" = "sex",
  "IMD" = "imd",
  "Ethnicity" = "ethnicity"
  )

readr::write_rds(
  list(demographic = demographic, clinical = clinical),
  here::here("output", "lib", "model_varlist.rds")
)

################################################################################
# strata vars for cox model ----
strata_vars <- c(
  "Region" = "region",
  "JCVI group" = "jcvi_group",
  "Date of eligibility for 1st dose" = "elig_date"
)

readr::write_rds(
  strata_vars,
  here::here("output", "lib", "strata_vars.rds")
)

################################################################################
# outcomes ----
outcomes <- c(
  "Any SARS-CoV-2 test" = "anytest", 
  "Positive SARS-CoV-2 test" = "postest", 
  "COVID-19 hospital admission" = "covidadmitted",
  "COVID-19 death" = "coviddeath", 
  "Non-COVID-19 death" = "noncoviddeath")

readr::write_rds(
  outcomes,
  here::here("output", "lib", "outcomes.rds")
)

################################################################################
# subgroups ----
subgroups <- c("16-64 years and clinically vulnerable", "18-39 years", "40-64 years", "65+ years")

readr::write_rds(
  subgroups,
  here::here("output", "lib", "subgroups.rds")
)
