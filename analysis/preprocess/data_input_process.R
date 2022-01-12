################################################################################

# This script:
# - reads, processes the extracted data, saves the following:
# - data_covs.rds = wide covariates and outcome data
# - data_*_vax_dates.rds = long and wide vaccine data
# - data_long_*_dates.rds = long covariates and outcome data

################################################################################

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

################################################################################
# inital pre-processing
cat("#### extract data ####\n")
data_extract <- 
  arrow::read_feather(file = here::here("output", "input.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(across(c(contains("_date"), dob), 
                ~ floor_date(
                  as.Date(., format="%Y-%m-%d"),
                  unit = "days"))) %>%
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
    ),
    #Subgroup
    subgroup = fct_case_when(
      jcvi_group %in% c("04", "06") & age_1 < 65 ~ "16-64 and clinically vulnerable",
      18 <= age_1 & age_1 < 40 ~ "18-39",
      40 <= age_1 & age_1 < 65 ~ "40-64",
      65 <= age_1 ~ "65+",
      TRUE ~ NA_character_
    )
    
  ) %>%
  select(-ethnicity_6, -ethnicity_6_sus) %>%
  droplevels()

################################################################################
cat("#### properties of data_processed ####\n")
# for checking for errors
data_properties(
  data = data_processed,
  path = file.path("output", "tables")
)

cat("## check subgroups as desired ##\n")
data_processed %>%
  group_by(subgroup, jcvi_group) %>%
  count() %>%
  ungroup()

################################################################################
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

# save long and wide datasets or vaccine variables
readr::write_rds(
  data_vax, 
  here::here("output", "data", "data_long_vax_dates.rds"), 
  compress="gz")

readr::write_rds(
  data_vax_wide,
  here::here("output", "data", "data_wide_vax_dates.rds"), 
  compress="gz")

###############################################################################
## create one-row-per-event datasets for recurring variables

# shielded
data_pr_shielded <- data_processed %>%
  select(patient_id,
         matches("^shielded\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "shielded_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

readr::write_rds(
  data_pr_shielded, 
  here::here("output", "data", "data_long_shielded_dates.rds"), 
  compress="gz")

###############################################################################
# nonshielded
data_pr_nonshielded <- data_processed %>%
  select(patient_id,
         matches("^nonshielded\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "nonshielded_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

readr::write_rds(
  data_pr_nonshielded, 
  here::here("output", "data", "data_long_nonshielded_dates.rds"), 
  compress="gz")

###############################################################################
# bmi
data_pr_bmi <- data_processed %>%
  select(patient_id,
         matches("^bmi\\_\\d+")) %>%
  rename_at(vars(contains("date")),
            ~ str_c("date_", str_extract(.x, "\\d+"))) %>%
  pivot_longer(
    cols = -patient_id,
    names_sep = "_",
    names_to = c(".value", "bmi_index"),
    values_drop_na = TRUE) %>%
  mutate(bmi = fct_case_when(
    bmi < 30 | bmi >=100 ~ "Not obese", # this cat includes missing and clinically implausible values
    bmi >= 30 & bmi < 35 ~ "Obese I (30-34.9)",
    bmi >= 35 & bmi < 40 ~ "Obese II (35-39.9)",
    bmi >= 40 & bmi < 100 ~ "Obese III (40+)",
    TRUE ~ NA_character_
  )) %>%
  arrange(patient_id, date) 

readr::write_rds(
  data_pr_bmi, 
  here::here("output", "data", "data_long_bmi_dates.rds"), 
  compress="gz")

###############################################################################
# suspected covid
data_pr_suspected_covid <- data_processed %>%
  select(patient_id,
         matches("^primary\\_care\\_suspected\\_covid\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "suspected_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

readr::write_rds(
  data_pr_suspected_covid, 
  here::here("output", "data", "data_long_pr_suspected_covid_dates.rds"), 
  compress="gz")

###############################################################################
# probable covid
data_pr_probable_covid <- data_processed %>%
  select(patient_id,
         matches("^primary\\_care\\_covid\\_case\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "probable_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

readr::write_rds(
  data_pr_probable_covid, 
  here::here("output", "data", "data_long_pr_probable_covid_dates.rds"),
  compress="gz")

###############################################################################
# covid admission
data_covidadmitted <- data_processed %>%
  select(patient_id, 
         matches("^covidadmitted\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "covidadmitted_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

readr::write_rds(
  data_covidadmitted, 
  here::here("output", "data", "data_long_covidadmitted_dates.rds"), 
  compress="gz")

###############################################################################
# positive test
data_postest <- data_processed %>%
  select(patient_id, 
         matches("^positive\\_test\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "postest_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

# combine outcomes where postest missing
# individuals with covidadmitted but not postest
data_covidadmitted_impute <- data_covidadmitted %>% 
  anti_join(data_postest, by = "patient_id") %>%
  rename(postest_index = covidadmitted_index)

# individuals with coviddeath but not postest
data_coviddeath_impute <- data_processed %>%
  select(patient_id, coviddeath_date) %>%
  filter(!is.na(coviddeath_date)) %>%
  anti_join(data_postest,
            by = "patient_id") %>%
  anti_join(data_covidadmitted_impute,
            by = "patient_id") %>%
  mutate(postest_index = "0") %>%
  rename(date = coviddeath_date)

#### may have to re-think this approach to combining outcomes if outcomes do 
#### become recurring rather than "ever"
data_postest <- bind_rows(
  data_postest,
  data_covidadmitted_impute,
  data_coviddeath_impute
)

readr::write_rds(
  data_postest, 
  here::here("output", "data", "data_long_postest_dates.rds"), 
  compress="gz")

################################################################################
# create dataset which contains the earliest date of any evidence of covid
# (not including covid death, as only applied to alive individuals)
data_covid_any <- list(
  data_pr_suspected_covid,
  data_pr_probable_covid,
  data_covidadmitted,
  data_postest
)

data_covid_any <- bind_rows(
  lapply(
    data_covid_any, 
    function(x) {
      name <- str_remove(names(x %>% select(ends_with("index"))), "_index")
      x %>% select(patient_id, date) %>% mutate(covid_event = name)
    }
    )
  ) %>%
  mutate(across(covid_event,
                factor,
                # if multiple recorded on the same date, this is the order of preference
                levels = c(
                  "covidadmitted",
                  "postest",
                  "probable",
                  "suspected"
                ))) %>%
  arrange(patient_id, date, covid_event) %>%
  # keep the first event to occur
  distinct(patient_id, .keep_all = TRUE) %>%
  rename(covid_any_date = date)

readr::write_rds(
  data_covid_any, 
  here::here("output", "data", "data_covid_any.rds"), 
  compress="gz")

################################################################################
# save dataset of covariates 
# (remove variables that are saved elsewhere)
readr::write_rds(
  data_processed %>%
    # remove vaccine variables
    select(-contains("_vax_")) %>%
    # remove recurring variables
    select(-starts_with(c(
      "shielded",
      "nonshielded",
      "bmi", 
      "primary_care_suspected_covid", 
      "primary_care_covid_case",
      "positive_test",
      "covidadmitted"))),
  here::here("output", "data", "data_covs.rds"), 
  compress="gz")
