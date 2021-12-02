library(tidyverse)
library(lubridate)
library(glue)

study_parameters <- readr::read_rds(here::here("output", "lib", "study_parameters.rds"))

set.seed(study_parameters$seed)

n <- study_parameters$n

start_date <- study_parameters$start_date

# input_old <- arrow::read_feather(file = here::here("output", "input_vax.feather")) 

regions <- readr::read_csv(here::here("output", "lib", "regions.csv"))

dummy_data <- tibble(
  
  patient_id = 1:n,
  
  age_1 = as.integer(runif(n, 16, 90), 0),
  
  sex = factor(rbernoulli(n, p=0.49),
               levels = c(FALSE, TRUE),
               labels = c("F", "M")),
  
  bmi_0 = rnorm(n, 25, 5),
  
  ethnicity_6 = factor(sample(
    x = c(as.character(1:5), NA_character_),
    size = n,
    replace = TRUE,
    prob = c(rep(0.99/5,5), 0.01)
  )),
  ethnicity_6_sus = factor(sample(
    x = c(as.character(1:5), NA_character_),
    size = n,
    replace = TRUE,
    prob = c(rep(0.99/5,5), 0.01)
  )),
  
  endoflife_0 = rbernoulli(n, p=0.005),
  midazolam_0 = rbernoulli(n, p=0.005),
  positive_test_0 = rbernoulli(n, p=0.05),
  primary_care_covid_case_0 = rbernoulli(n, p=0.05),
  primary_care_suspected_covid_0 = rbernoulli(n, p=0.05),
  covidadmitted_0 = rbernoulli(n, p=0.05),
  
  endoflife_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  midazolam_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  positive_test_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  primary_care_covid_case_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  primary_care_suspected_covid_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  covidadmitted_0_date = as_date(start_date) + days(sample(x = 1:600, size = n, replace = TRUE)),
  
  hscworker = rbernoulli(n, p=0.01),
  
  imd_0 = sample(
    x = seq(100L,32100L,100L),
    size = n,
    replace = TRUE
  ),
  
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
  mutate(covid_vax_pfizer_1_date = elig_date + days(round(rnorm(nrow(.), mean = 10, sd = 3))),
         covid_vax_az_1_date = elig_date + days(round(rnorm(nrow(.), mean = 10, sd = 3))),
         covid_vax_moderna_1_date = elig_date + days(round(rnorm(nrow(.), mean = 10, sd = 3)))) %>%
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

arrow::write_feather(dummy_data, here::here("analysis", "lib", "dummy_data_vax.feather"))

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
