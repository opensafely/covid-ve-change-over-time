################################################################################

# This script:
# - reads the processed data
# - applies eligibility criteria from boxes a and b of Figure 3 in protocol

################################################################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

# read processed covariates data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

# read wide vaccine dates data
data_vax_wide <- readr::read_rds(
  here::here("output", "data", "data_wide_vax_dates.rds"))

# count the number of patients in the extracted data
eligibility_count <- tribble(
  ~description, ~n,
  "Extracted using study_definition", n_distinct(data_processed$patient_id)
)

################################################################################
cat("#### apply exclusion criteria from box a to processed data ####\n")
# remove dummy jcvi_group
data_eligible_a <- data_processed %>%
  filter(!(jcvi_group %in% "99"))

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Samples with JCVI group 99 removed",
    n =  n_distinct(data_eligible_a$patient_id)
  )

# remove jcvi_group=01 (care homes)
data_eligible_a <- data_eligible_a %>%
  filter(!(jcvi_group %in% "01"))

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Samples with JCVI group 01 removed",
    n =  n_distinct(data_eligible_a$patient_id)
  )

# remove if aged under 16 or over 120 (the latter to avoid probable errors)
data_eligible_a <- data_eligible_a %>%
  filter(age_1 >= 16, age_1 <= 120)

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Samples with age_1 < 16 and > 120 removed.",
    n =  n_distinct(data_eligible_a$patient_id)
  )

# remove if any missing data for key variables
data_eligible_a <- data_eligible_a %>%
  filter(
    !is.na(sex),
    !is.na(region_0),
    !is.na(ethnicity),
    !is.na(imd_0))

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Samples with missing ethnicity, sex, imd, region removed.",
    n =  n_distinct(data_eligible_a$patient_id)
  )

# remove if evidence of covid infection on or before elig_date + 42 days
# COVID admission
data_eligible_a <- data_eligible_a %>%
  filter(
    ! (
      !is.na(covid_any_date) &
        (covid_any_date <= elig_date + days(42)) &
        covid_event %in% "covidadmitted"
    ))

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Samples with prior COVID admission removed.",
    n =  n_distinct(data_eligible_a$patient_id)
  )

# positive COVID test
data_eligible_a <- data_eligible_a %>%
  filter(
    ! (
      !is.na(covid_any_date) &
        (covid_any_date <= elig_date + days(42)) &
        covid_event %in% "postest"
    ))

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Samples with prior positive COVID test removed.",
    n =  n_distinct(data_eligible_a$patient_id)
  )

# probable COVID 
data_eligible_a <- data_eligible_a %>%
  filter(
    ! (
      !is.na(covid_any_date) &
        (covid_any_date <= elig_date + days(42)) &
        covid_event %in% "probable"
    ))

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Samples with prior probable COVID removed.",
    n =  n_distinct(data_eligible_a$patient_id)
  )

# suspected COVID 
data_eligible_a <- data_eligible_a %>%
  filter(
    ! (
      !is.na(covid_any_date) &
        (covid_any_date <= elig_date + days(42)) &
        covid_event %in% "suspected"
    ))

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Samples with prior suspected COVID removed.",
    n =  n_distinct(data_eligible_a$patient_id)
  )

readr::write_rds(data_eligible_a %>%
                   select(patient_id, jcvi_group, elig_date, region_0, ethnicity) %>%
                   droplevels(),
                 here::here("output", "data", "data_eligible_a.rds"),
                 compress="gz")

################################################################################
cat("#### apply exclusion criteria from box b to processed data ####\n")
data_eligible_b <- data_eligible_a %>%
  left_join(data_vax_wide, by = "patient_id") %>%
  mutate(
    between_doses = as.numeric(covid_vax_2_date - covid_vax_1_date)/7
  ) %>%
  filter(
    ### inclusion
    # first and second dose same brand
    covid_vax_1_brand == covid_vax_2_brand,
    # second dose must be az or pfizer
    covid_vax_1_brand %in% c("az", "pfizer"),
    # second dose received in interval [elig_date + 6 weeks, elig_date + 16 weeks)
    covid_vax_2_date >= elig_date + weeks(6),
    covid_vax_2_date < elig_date + weeks(16),
    
    ### exclusion
    # first dose received before eligibility date
    covid_vax_1_date > elig_date,
    # less than six or more than 14 weeks between 1st and 2nd dose
    between_doses >= 6,
    between_doses < 14,
    # flagged as hcw
    !hscworker
  ) %>%
  select(patient_id, jcvi_group, elig_date, region_0, ethnicity) %>%
  droplevels()

eligibility_count <- eligibility_count %>%
  add_row(
    description = "Criteria in box b applied.",
    n = n_distinct(data_eligible_b$patient_id)
  )

readr::write_rds(data_eligible_b,
                 here::here("output", "data", "data_eligible_b.rds"),
                 compress="gz")

# number of people eligible at each stage ----
eligibility_count <- eligibility_count %>%
  # round to nearest 10
  mutate(across(n, ~round(.x, -1))) %>%
  mutate(n_removed = lag(n) - n)

readr::write_csv(eligibility_count,
                 here::here("output", "lib", "eligibility_count_ab.csv"))


################################################################################
# jcvi_group, elig_date combos ----
fix_age <- data_processed %>%
  mutate(age = if_else(
    jcvi_group %in% c("10", "11", "12"),
    age_2,
    age_1)) %>%
  select(patient_id, age)

group_age_ranges <- data_eligible_b %>%
  left_join(fix_age, 
            by = "patient_id") %>%
  group_by(jcvi_group, elig_date) %>%
  summarise(min = min(age), max = max(age), .groups = "keep") %>%
  ungroup() 

# check none of the min / max correspond to < 5 individuals
check <- function(x) {
  fix_age %>%
    filter(age %in% x) %>%
    distinct(patient_id) %>%
    nrow(.)
}

group_age_ranges <- group_age_ranges %>%
  mutate(
    check_min = sapply(
      group_age_ranges$min,
      check),
    check_max = sapply(
      group_age_ranges$max,
      check)
  ) %>%
  mutate(age_range = case_when(
    max > 80 ~ glue("{min} +"),
    check_min < 5 | check_max < 5 ~ "[REDACTED]",
    TRUE ~ glue("{min} - {max}")
  )) %>%
  mutate(across(age_range, as.character)) %>%
  select(-ends_with(c("min", "max")))

# save as csv so that it can be checked
readr::write_csv(group_age_ranges,
                 here::here("output", "lib", "group_age_ranges.csv"))