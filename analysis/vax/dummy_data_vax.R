######################################

# This script:
# - creates dummy data for study_definition_vax.py

######################################

library(tidyverse)
library(lubridate)
library(glue)

source(here::here("analysis", "lib", "dummy_data_functions.R"))

vars_date_0 <- c("endoflife_0_date",
                 "midazolam_0_date", 
                 "positive_test_0_date", 
                 "primary_care_covid_case_0_date", 
                 "primary_care_suspected_covid_0_date", 
                 "covidadmitted_0_date")


jcvi_group_patterns <- readr::read_csv(here::here("output", "lib", "jcvi_groups.csv")) %>%
  mutate(across(definition, ~str_extract(.x, "age_. >=\\d{2}"))) %>%
  mutate(across(definition, ~case_when(group=="01" ~ "age_1 >=90",
                                       group=="06" ~ "age_1 >=62",
                                       !is.na(.x) ~ .x,
                                       TRUE ~ "TRUE")))

# conditions for eligibility dates
elig_date_patterns <- readr::read_csv(here::here("output", "lib", "elig_dates.csv")) %>%
  mutate(across(description, ~str_replace_all(.x, "p=", "p=="))) %>%
  mutate(across(description, ~str_replace_all(.x, "OR", "|"))) %>%
  mutate(across(description, ~str_replace_all(.x, "AND", "&"))) %>%
  mutate(across(description, ~str_replace_all(.x, "DEFAULT", "TRUE")))

dummy_data <- tibble(patient_id = 1:n) %>%
  mutate(age_1 = as.integer(runif(nrow(.), 16, 90), 0),
         age_2 = age_1,
         bmi_0 = rnorm(n, 25, 5),
         imd_0 = sample(
           x = seq(100L,32100L,100L),
           size = n,
           replace = TRUE)) %>%
  var_category(sex, categories = c("F", "M")) %>%
  var_category(ethnicity_6, 
               categories = c(as.character(1:5), NA_character_),
               ratios = c(rep(0.99/5,5), 0.01)) %>%
  var_category(ethnicity_6_sus, 
               categories = c(as.character(1:5), NA_character_),
               ratios = c(rep(0.99/5,5), 0.01)) %>%
  var_category(region_0, categories = regions$region, ratios = regions$ratio) %>%
  # binary vars for exclusion criteria
  bind_cols(
    pmap(
      list(a = c("hscworker"), 
           b = c(0.01)),
      function(a,b) 
        var_binary(
          .data = ., 
          name = !! a,
          incidence = b,
          keep_vars = FALSE
        ))) %>%
  # date vars for exclusion criteria
  bind_cols(
    pmap(
      list(a = vars_date_0, 
           b = rep(0.01, length(vars_date_0))),
      function(a,b) 
        var_date(
        .data = ., 
        name = !! a,
        incidence = b,
        keep_vars = FALSE
      ))) %>%
  # jcvi_group
  var_category(
    name = jcvi_group, 
    categories = jcvi_group_patterns$group, 
    conditions = jcvi_group_patterns$definition
    ) %>%
  # elig_date
  var_category(
    name = elig_date,
    categories = elig_date_patterns$date,
    conditions = elig_date_patterns$description)
  

# fix vaccine dates so that they have roughly correct distribution
dummy_data <- dummy_data %>%
  mutate(
    covid_vax_pfizer_1_date = as.Date(elig_date) + days(round(rnorm(nrow(.), mean = 10, sd = 3))),
    covid_vax_az_1_date = as.Date(elig_date) + days(round(rnorm(nrow(.), mean = 10, sd = 3))),
    covid_vax_moderna_1_date = as.Date(elig_date) + days(round(rnorm(nrow(.), mean = 10, sd = 3)))) %>%
  mutate(
    vaccine_1_type = sample(
      x = c("pfizer", "az", "moderna", "none"), 
      size = nrow(.),
      replace = TRUE,
      prob = c(0.4, 0.4, 0.1, 0.1)
      ),
    missing_pfizer_2 = rbernoulli(nrow(.), p=0.05),
    missing_az_2 = rbernoulli(nrow(.), p=0.05),
    missing_moderna_2 = rbernoulli(nrow(.), p=0.05),
    missing_pfizer_3 = rbernoulli(nrow(.), p=0.9),
    missing_az_3 = rbernoulli(nrow(.), p=0.9),
    missing_moderna_3 = rbernoulli(nrow(.), p=0.9)
  ) %>%
  mutate(across(covid_vax_pfizer_1_date, 
                ~if_else(
                  vaccine_1_type %in% "pfizer",
                  .x,
                  NA_Date_))) %>%
  mutate(across(covid_vax_az_1_date, 
                ~if_else(
                  vaccine_1_type %in% "az",
                  .x,
                  NA_Date_))) %>%
  mutate(across(covid_vax_moderna_1_date, 
                ~if_else(
                  vaccine_1_type %in% "moderna",
                  .x,
                  NA_Date_))) %>%
  mutate(across(matches("covid_vax_\\w+_1_date"),
                ~ if_else(
                  vaccine_1_type %in% "none",
                  NA_Date_,
                  .x
                ))) %>%
  mutate(
    covid_vax_pfizer_2_date = covid_vax_pfizer_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 3))),
    covid_vax_az_2_date = covid_vax_az_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 3))),
    covid_vax_moderna_2_date = covid_vax_moderna_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 3))),
    ) %>%
  mutate(across(covid_vax_pfizer_2_date, 
                ~if_else(
                  missing_pfizer_2,
                  NA_Date_,
                  .x))) %>%
  mutate(across(covid_vax_az_2_date, 
                ~if_else(
                  missing_az_2,
                  NA_Date_,
                  .x))) %>%
  mutate(across(covid_vax_moderna_2_date, 
                ~if_else(
                  missing_moderna_2,
                  NA_Date_,
                  .x))) %>%
  mutate(
    covid_vax_pfizer_3_date = covid_vax_pfizer_2_date + days(round(rnorm(nrow(.), mean = 6*4*7, sd = 7))),
    covid_vax_az_3_date = covid_vax_az_2_date + days(round(rnorm(nrow(.), mean = 6*4*7, sd = 7))),
    covid_vax_moderna_3_date = covid_vax_moderna_2_date + days(round(rnorm(nrow(.), mean = 6*4*7, sd = 7))),
  ) %>%
  select(-starts_with("missing"), -vaccine_1_type) %>%
  mutate(across(ends_with("date"), as.POSIXct))

arrow::write_feather(dummy_data, here::here("analysis", "vax", "dummy_data_vax.feather"))

# #  checks
# # all names there and the same?
# all(sort(names(dummy_data)) == sort(names(input_old)))
# # any different types
# classes_input_old <- sapply(
#   sort(names(input_old)),
#   function(x) 
#     class(input_old[[x]]))
# classes_dummy_data <- sapply(
#   sort(names(dummy_data)), 
#   function(x) 
#     class(dummy_data[[x]]))
# classes_match <- sapply(
#   names(classes_input_old),
#   function(x)
#     all(classes_input_old[[x]] == classes_dummy_data[[x]]))

# sort(names(dummy_data))[!classes_match]
