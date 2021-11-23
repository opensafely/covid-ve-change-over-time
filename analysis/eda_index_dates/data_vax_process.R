######################################

# This script:
# - reads the extracted data
# - processes the extracted data
# - cleans the vaccination data to identify second doses of pfizer and az
# - creates a dataset with the (redacted) number of individuals receiving their second vaccination on each date in a sequence
# - (the date sequence depends on their vaccine eligibility date and counts are stratified by region and vaccine brand)

######################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

## source functions
source(here::here("analysis", "lib", "data_properties.R"))

## create folders for outputs
data_dir <- here::here("output", "data")
eda_data_dir <- here::here("output", "eda_index_dates", "data")
dir.create(data_dir, showWarnings = FALSE, recursive=TRUE)
dir.create(eda_data_dir, showWarnings = FALSE, recursive=TRUE)

## import dates
dates <- readr::read_rds(here::here("output", "lib", "study_dates.rds"))

# Custom functions
fct_case_when <- function(...) {
  # uses dplyr::case_when but converts the output to a factor,
  # with factors ordered as they appear in the case_when's  ... argument
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}

cat("#### print variable names ####\n")
arrow::read_feather(here::here("output", "input_vax.feather")) %>%
  names() %>%
  sort() %>%
  print()

cat("#### extract data ####\n")
data_extract <- 
  arrow::read_feather(file = here::here("output", "input_vax.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(across(contains("_date"), ~ as.Date(., format="%Y-%m-%d"))) 

cat("#### check format of elig_date ####\n")
elig_date_test <- data_extract %>%
  distinct(elig_date)

# only for dummy data:
if (nrow(elig_date_test) <= 2) {
  
  cat("#### fix dummy data ####\n")
  jcvi_groups <- readr::read_csv(here::here("output", "lib", "jcvi_groups.csv"))
  
  elig_dates <- readr::read_csv(here::here("output", "lib", "elig_dates.csv"))
  
  jcvi_group_cases <- jcvi_groups %>%
    mutate(across(definition, ~str_extract(.x, "age_. >=\\d{2}"))) %>%
    mutate(across(definition, ~case_when(group=="01" ~ "age_1 >=90",
                                         group=="06" ~ "age_1 >=62",
                                         !is.na(.x) ~ .x,
                                         TRUE ~ "TRUE"))) %>%
    transmute(cases = str_c(definition, " ~ \'", group, "\'")) %>%
    unlist() %>% 
    unname() %>%
    str_c(., collapse = ", ")
  
  elig_date_cases <- elig_dates %>%
    mutate(across(description, ~str_replace_all(.x, "p=", "p=="))) %>%
    mutate(across(description, ~str_replace_all(.x, "OR", "|"))) %>%
    mutate(across(description, ~str_replace_all(.x, "AND", "&"))) %>%
    mutate(across(description, ~str_replace_all(.x, "DEFAULT", "TRUE"))) %>%
    transmute(cases = str_c(description, " ~ \'", date, "\'")) %>%
    unlist() %>% 
    unname() %>%
    str_c(., collapse = ", ")
  
  # JCVI groups based on age
  data_extract <- eval(parse(text = glue("data_extract %>% mutate(jcvi_group = case_when({jcvi_group_cases}))")))
  # eligibility dates based on JCVI groups
  data_extract <- eval(parse(text = glue("data_extract %>% mutate(elig_date = as.Date(case_when({elig_date_cases}), format=\'%Y-%m-%d\'))")))
  
  # fix vaccine dates so that they have roughly correct distribution
  data_extract <- data_extract %>%
    mutate(r1 = round(rnorm(nrow(.), mean = 7, sd = 7)),
           r2 = round(rnorm(nrow(.), mean = 7, sd = 7)),
           r3 = round(rnorm(nrow(.), mean = 10*7, sd = 7)),
           r4 = round(rnorm(nrow(.), mean = 10*7, sd = 7)),
           r5 = rbernoulli(nrow(.), p=0.01),
           r6 = rbernoulli(nrow(.), p=0.01)) %>%
    mutate(across(covid_vax_pfizer_1_date, ~if_else(!is.na(.x), elig_date + days(r1),.x))) %>%
    mutate(across(covid_vax_az_1_date, ~if_else(!is.na(.x), elig_date + days(r2),.x))) %>%
    mutate(covid_vax_pfizer_2_date = covid_vax_pfizer_1_date + days(r3),
           covid_vax_az_2_date = covid_vax_az_1_date + days(r4)) %>%
    mutate(across(contains("vax_moderna"), ~if_else(r5, .x, NA_Date_))) %>%
    mutate(across(contains("vax_disease"), ~if_else(r6, .x, NA_Date_))) %>%
    select(-matches("r\\d{1}"))
  
}

cat("#### process extracted data ####\n")
data_vax_processed <- data_extract %>%
  # derive ethnicity variable
  mutate(
    # Ethnicity
    ethnicity = if_else(is.na(ethnicity_6), ethnicity_6_sus, ethnicity_6),
    ethnicity = fct_case_when(
      ethnicity == "1" ~ "White",
      ethnicity == "4" ~ "Black",
      ethnicity == "3" ~ "South Asian",
      ethnicity == "2" ~ "Mixed",
      ethnicity == "5" ~ "Other",
      TRUE ~ NA_character_
    ),
    # IMD quintile
    imd_0 = fct_case_when(
      imd_0 == 1 ~ "1 most deprived",
      imd_0 == 2 ~ "2",
      imd_0 == 3 ~ "3",
      imd_0 == 4 ~ "4",
      imd_0 == 5 ~ "5 least deprived",
      TRUE ~ NA_character_
    ),
    # Sex
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      #sex == "I" ~ "Inter-sex",
      #sex == "U" ~ "Unknown",
      TRUE ~ NA_character_
    ),
    # BMI
    bmi_0 = fct_case_when(
      bmi_0 < 30 | bmi_0 >=100 ~ "Not obese", # this cat includes missing and clinically implausible values
      bmi_0 >= 30 & bmi_0 < 35 ~ "Obese I (30-34.9)",
      bmi_0 >= 35 & bmi_0 < 40 ~ "Obese II (35-39.9)",
      bmi_0 >= 40 & bmi_0 < 100 ~ "Obese III (40+)",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-ethnicity_6, -ethnicity_6_sus) %>%
  droplevels()

cat("#### properties of data_vax_processed ####\n")
data_properties(
  data = data_vax_processed,
  path = data_dir
)  

cat("#### apply exclusion criteria to processed data ####\n")
data_eligible <- data_vax_processed %>%
  # apply exclusion criteria
  filter(
    # remove if any missing data for key variables
    !is.na(ethnicity),
    !is.na(sex),
    !is.na(imd_0),
    !is.na(region_0),
    # remove if initiated end of life care on or before elig_date + 42 days
    is.na(endoflife_0_date),
    is.na(midazolam_0_date),
    # remove if evidence of covid infection on or before elig_date + 42 days
    is.na(positive_test_0_date),
    is.na(primary_care_covid_case_0_date),
    is.na(primary_care_suspected_covid_0_date),
    is.na(covidadmitted_0_date),
    # remove dummy groups and dates
    !(jcvi_group %in% "99"),
    !(elig_date %in% as.Date("2100-12-31"))
  ) 

# process vaccine data
data_vax <- local({
  
  data_vax_pfizer <- data_eligible %>%
    select(patient_id, matches("covid\\_vax\\_pfizer\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_pfizer_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_az <- data_eligible %>%
    select(patient_id, matches("covid\\_vax\\_az\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_az_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_moderna <- data_eligible %>%
    select(patient_id, matches("covid\\_vax\\_moderna\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_moderna_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  
  data_vax <-
    data_vax_pfizer %>%
    full_join(data_vax_az, by=c("patient_id", "date")) %>%
    full_join(data_vax_moderna, by=c("patient_id", "date")) %>%
    mutate(
      brand = fct_case_when(
        (!is.na(vax_az_index)) & is.na(vax_pfizer_index) & is.na(vax_moderna_index) ~ "az",
        is.na(vax_az_index) & (!is.na(vax_pfizer_index)) & is.na(vax_moderna_index) ~ "pfizer",
        is.na(vax_az_index) & is.na(vax_pfizer_index) & (!is.na(vax_moderna_index)) ~ "moderna",
        (!is.na(vax_az_index)) + (!is.na(vax_pfizer_index)) + (!is.na(vax_moderna_index)) > 1 ~ "duplicate",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(patient_id, date) %>%
    group_by(patient_id) %>%
    mutate(
      vax_index=row_number()
    ) %>%
    ungroup()
  
  data_vax
  
})

data_vax_wide <- data_vax %>%
  pivot_wider(
    id_cols= patient_id,
    names_from = c("vax_index"),
    values_from = c("date", "brand"),
    names_glue = "covid_vax_{vax_index}_{.value}"
  )

readr::write_rds(data_vax, 
                 here::here("output", "data", "data_long_vax_dates.rds"), 
                 compress="gz")
readr::write_rds(data_vax_wide,
          here::here("output", "data", "data_wide_vax_dates.rds"), 
          compress="gz")

cat("#### derive data_2nd_dose ####\n")
data_2nd_dose <- data_vax %>%
  # keep first dose if first dose of pfizer or az
  # keep second dose if second dose of pfizer or az
  filter(
    vax_index %in% c(1,2),
    vax_index == vax_pfizer_index | vax_index == vax_az_index
  ) %>%
  select(patient_id, vax_index, brand, date) %>%
  pivot_wider(names_from = vax_index,
              values_from = date, 
              names_prefix = "dose_") %>%
  left_join(
    data_eligible %>%
      select(patient_id, jcvi_group, elig_date, region_0, hscworker),
    by = "patient_id") %>%
  # first dose must have occurred on or after elig_date
  filter(dose_1 >= elig_date) %>%
  # only keep if 2nd dose received [6,14) weeks after 1st dose
  filter(dose_1 + weeks(6) <= dose_2 & dose_2 < dose_1 + weeks(14)) %>%
  filter(!hscworker) %>%
  select(-hscworker)
  
  
readr::write_rds(data_2nd_dose,
                 here::here("output", "eda_index_dates", "data", "data_2nd_dose.rds"), 
                 compress="gz")

