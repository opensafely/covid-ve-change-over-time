######################################

# This script:
# - reads the extracted data
# - processes the extracted data
# - saves processed data

######################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

## source functions
source(here::here("analysis", "lib", "data_process_functions.R"))
source(here::here("analysis", "lib", "data_properties.R"))

## create folders for outputs
dir.create(here::here("output", "data"), showWarnings = FALSE, recursive=TRUE)
dir.create(here::here("output", "tables"), showWarnings = FALSE, recursive=TRUE)

## import study_parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

## regions
regions <- readr::read_csv(
  here::here("output", "lib", "regions.csv")
)

cat("#### extract data ####\n")
data_extract <- 
  arrow::read_feather(file = here::here("output", "input.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(across(contains("_date"), ~ as.Date(., format="%Y-%m-%d"))) %>%
  mutate(across(imd_0, ~as.integer(as.character(.x))))

cat("#### process extracted data ####\n")
data_processed <- data_extract %>%
  # derive ethnicity variable
  mutate(
    # Region
    region_0 = factor(region_0, levels = regions$region),
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
      imd_0 < 1 | is.na(imd_0) ~ NA_character_,
      imd_0 < 32844*1/5 ~ "1 most deprived",
      imd_0 < 32844*2/5 ~ "2",
      imd_0 < 32844*3/5 ~ "3",
      imd_0 < 32844*4/5 ~ "4",
      TRUE ~ "5 least deprived"
    ),
    # Sex
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      #sex == "I" ~ "Inter-sex",
      #sex == "U" ~ "Unknown",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-ethnicity_6, -ethnicity_6_sus) %>%
  droplevels()

# process vaccine data
data_vax <- local({
  
  data_vax_pfizer <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_pfizer\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_pfizer_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_az <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_az\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_az_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_moderna <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_moderna\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_moderna_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  
  data_vax <- data_processed %>% # to get the unvaccinated
    # filter(if_all(starts_with("covid_vax"), ~ is.na(.))) %>%
    filter_at(vars(starts_with("covid_vax")), all_vars(is.na(.))) %>%
    select(patient_id) %>% 
    full_join(
      data_vax_pfizer %>%
        full_join(data_vax_az, by=c("patient_id", "date")) %>%
        full_join(data_vax_moderna, by=c("patient_id", "date")),
      by = "patient_id"
    ) %>%
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
    ungroup() %>%
    droplevels()
  
  data_vax
  
})

data_vax_wide <- data_vax %>%
  pivot_wider(
    id_cols= patient_id,
    names_from = c("vax_index"),
    values_from = c("date", "brand"),
    names_glue = "covid_vax_{vax_index}_{.value}"
  )


cat("#### properties of data_processed ####\n")
data_properties(
  data = data_processed,
  path = file.path("output", "tables")
)

# save dataset of covariates 
# (i.e. remove vaccine variables as they are saved elsewhere)
readr::write_rds(
  data_processed %>%
    select(-contains("_vax_")),
  here::here("output", "data", "data_covs.rds"), 
  compress="gz")

# save long and wide datasets or vaccine variables
readr::write_rds(
  data_vax, 
  here::here("output", "data", "data_long_vax_dates.rds"), 
  compress="gz")

readr::write_rds(
  data_vax_wide,
  here::here("output", "data", "data_wide_vax_dates.rds"), 
  compress="gz")
