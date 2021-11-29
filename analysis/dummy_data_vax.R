library(tidyverse)
library(lubridate)
library(glue)

set.seed(4543)

n <- 10000

start_date <- "2020-12-08"

input_old <- 
  arrow::read_feather(file = here::here("output", "input_vax.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(across(contains("_date"), ~ as.Date(., format="%Y-%m-%d"))) 

regions <- readr::read_csv(here::here("output", "lib", "regions.csv"))

dummy_data <- tibble(
  
  patient_id = 1:n,
  
  age_1 = as.integer(runif(n, 0,100),0),
  
  sex = factor(rbernoulli(n, p=0.49),
               levels = c(FALSE, TRUE),
               labels = c("F", "M")),
  
  bmi_0 = rnorm(n, 25, 5),
  
  ethnicity_6 = factor(sample(
    x = c(as.character(1:5), NA_character_),
    size = n,
    replace = TRUE
  )),
  ethnicity_6_sus = factor(sample(
    x = c(as.character(1:5), NA_character_),
    size = n,
    replace = TRUE
  )),
  
  endoflife_0 = rbernoulli(n, p=0.005),
  midazolam_0 = rbernoulli(n, p=0.005),
  positive_test_0 = rbernoulli(n, p=0.1),
  primary_care_covid_case_0 = rbernoulli(n, p=0.1),
  primary_care_suspected_covid_0 = rbernoulli(n, p=0.1),
  covidadmitted_0 = rbernoulli(n, p=0.1),
  
  endoflife_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  midazolam_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  positive_test_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  primary_care_covid_case_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  primary_care_suspected_covid_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  covidadmitted_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  
  hscworker = rbernoulli(n, p=0.01),
  
  imd_0 = factor(sample(
    x = as.character(seq(100,32100,100)),
    size = n,
    replace = TRUE
  )),
  
  region_0 = factor(sample(
    x = regions$region,
    size = n,
    replace = TRUE,
    prob = regions$ratio
  ))
  
) %>%
  mutate(
    age_2 = age_1,
    
    endoflife_0_date = if_else(endoflife_0, endoflife_0_date, NA_Date_),
    midazolam_0_date = if_else(midazolam_0, midazolam_0_date, NA_Date_),
    positive_test_0_date = if_else(positive_test_0, positive_test_0_date, NA_Date_),
    primary_care_covid_case_0_date = if_else(primary_care_covid_case_0, primary_care_covid_case_0_date, NA_Date_),
    primary_care_suspected_covid_0_date = if_else(primary_care_suspected_covid_0, primary_care_suspected_covid_0_date, NA_Date_),
    covidadmitted_0_date = if_else(covidadmitted_0, covidadmitted_0_date, NA_Date_)
    
  ) %>%
  select(-endoflife_0, -midazolam_0, -positive_test_0, -primary_care_covid_case_0, -primary_care_suspected_covid_0, -covidadmitted_0)


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
dummy_data <- eval(parse(text = glue("dummy_data %>% mutate(jcvi_group = factor(case_when({jcvi_group_cases})))")))
# eligibility dates based on JCVI groups
dummy_data <- eval(parse(text = glue("dummy_data %>% mutate(elig_date = as.Date(case_when({elig_date_cases}), format=\'%Y-%m-%d\'))")))


# fix vaccine dates so that they have roughly correct distribution
dummy_data <- dummy_data %>%
  mutate(covid_vax_pfizer_1_date = elig_date + days(round(rnorm(nrow(.), mean = 10, sd = 7))),
         covid_vax_az_1_date = elig_date + days(round(rnorm(nrow(.), mean = 10, sd = 7))),
         covid_vax_moderna_1_date = elig_date + days(round(rnorm(nrow(.), mean = 10, sd = 7)))) %>%
  mutate(
    missing_pfizer_1 = rbernoulli(nrow(.), p=0.3),
    missing_az_1 = rbernoulli(nrow(.), p=0.3),
    missing_moderna_1 = rbernoulli(nrow(.), p=0.9),
    missing_pfizer_2 = rbernoulli(nrow(.), p=0.3),
    missing_az_2 = rbernoulli(nrow(.), p=0.3),
    missing_moderna_2 = rbernoulli(nrow(.), p=0.9),
    missing_pfizer_3 = rbernoulli(nrow(.), p=0.3),
    missing_az_3 = rbernoulli(nrow(.), p=0.3),
    missing_moderna_3 = rbernoulli(nrow(.), p=0.9),
  ) %>%
  mutate(across(covid_vax_pfizer_1_date, 
                ~if_else(
                  missing_pfizer_1, 
                  NA_Date_,
                  .x))) %>%
  mutate(across(covid_vax_az_1_date, 
                ~if_else(
                  missing_az_1,
                  NA_Date_,
                  .x))) %>%
  mutate(across(covid_vax_moderna_1_date, 
                ~if_else(
                  missing_moderna_1,
                  NA_Date_,
                  .x))) %>%
  mutate(
    covid_vax_pfizer_2_date = covid_vax_pfizer_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 7))),
    covid_vax_az_2_date = covid_vax_az_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 7))),
    covid_vax_moderna_2_date = covid_vax_moderna_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 7))),
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
  select(-starts_with("missing"))

# readr::write_csv(dummy_data, here::here("analysis", "lib", "dummy_data.csv"))
arrow::write_feather(dummy_data, here::here("analysis", "lib", "dummy_data.feather"))

#  checks
# all names there and the same?
all(sort(names(dummy_data)) == sort(names(input_old)))
# any different types
sort(names(dummy_data))[sapply(
  sort(names(input_old)),
  function(x) 
    class(input_old[[x]])) != sapply(
      sort(names(dummy_data)), 
      function(x) 
        class(dummy_data[[x]]))]
