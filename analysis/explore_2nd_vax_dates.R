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
library(readr)
library(glue)

## source functions
source(here::here("analysis", "lib", "data_properties.R"))

## create folders for outputs
data_dir <- here::here("output", "explore_2nd_vax_dates", "data")
dir.create(data_dir, showWarnings = FALSE, recursive=TRUE)

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
read_csv(here::here("output", "input.csv"),
         n_max = 0,
         col_types = cols()) %>%
  names() %>%
  sort() %>%
  print()

cat("#### extract data ####\n")
data_extract0 <- read_csv(
  file = here::here("output", "input.csv"),
  col_types = cols_only(
    
    ## Identifier
    patient_id = col_integer(),
    
    ## Eligibility dates
    elig_date = col_date(format="%Y-%m-%d"),
    jcvi_group = col_character(),
    
    ## UK region
    region = col_character(),
    
    ## Age
    age_1 = col_integer(),
    age_2 = col_integer(),
    
    ## Variables for applying exclusion criteria
    ## COVID
    positive_test_0_date = col_date(format="%Y-%m-%d"),
    primary_care_covid_case_0_date = col_date(format="%Y-%m-%d"),
    primary_care_suspected_covid_0_date = col_date(format="%Y-%m-%d"),
    covidadmitted_0_date = col_date(format="%Y-%m-%d"),
    ## clinical
    longres_0_date = col_date(format="%Y-%m-%d"),
    endoflife_0_date = col_date(format="%Y-%m-%d"),
    midazolam_0_date = col_date(format="%Y-%m-%d"),
    ## missing
    sex = col_character(),
    ethnicity_6 = col_character(),
    ethnicity_6_sus = col_character(),
    imd = col_character(),
    
    ## Vaccination variables
    # First COVID vaccination date
    covid_vax_az_1_date = col_date(format="%Y-%m-%d"),
    covid_vax_az_2_date = col_date(format="%Y-%m-%d"),
    covid_vax_az_3_date = col_date(format="%Y-%m-%d"),
    covid_vax_disease_1_date = col_date(format="%Y-%m-%d"),
    covid_vax_disease_2_date = col_date(format="%Y-%m-%d"),
    covid_vax_disease_3_date = col_date(format="%Y-%m-%d"),
    covid_vax_moderna_1_date = col_date(format="%Y-%m-%d"),
    covid_vax_moderna_2_date = col_date(format="%Y-%m-%d"),
    covid_vax_moderna_3_date = col_date(format="%Y-%m-%d"),
    covid_vax_pfizer_1_date = col_date(format="%Y-%m-%d"),
    covid_vax_pfizer_2_date = col_date(format="%Y-%m-%d"),
    covid_vax_pfizer_3_date = col_date(format="%Y-%m-%d")
    
  ),
  na = character() # more stable to convert to missing later
) 

cat("#### parse NAs ####\n")
data_extract <- data_extract0 %>%
  mutate(across(
    .cols = where(is.character),
    .fns = ~na_if(.x, "")
  )) %>%
  mutate(patient_id = row_number()) %>% # create new ID variable, as duplicates after binding
  arrange(patient_id)

cat("#### check format of elig_date ####\n")
elig_date_test <- data_extract %>%
  select(elig_date) %>%
  filter(!is.na(elig_date) &
           str_detect(as.character(elig_date), "\\d{4}-\\d{2}-\\d{2}"))

# only for dummy data:
if (nrow(elig_date_test) == 0) {
  
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
data_processed <- data_extract %>%
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
    # imd quintile
    imd = na_if(imd, "0"),
    imd = fct_case_when(
      imd == 1 ~ "1 most deprived",
      imd == 2 ~ "2",
      imd == 3 ~ "3",
      imd == 4 ~ "4",
      imd == 5 ~ "5 least deprived",
      TRUE ~ NA_character_
    )
    
  ) %>%
  select(-ethnicity_6, -ethnicity_6_sus)

cat("#### properties of data_processed ####\n")
data_properties(
  data = data_processed,
  path = data_dir
)  

cat("#### apply exclusion criteria to processed data ####\n")
data_eligible <- data_processed %>%
  # apply exclusion criteria
  filter(
    # remove if any missing data for key variables
    !is.na(ethnicity),
    !is.na(sex),
    !is.na(imd),
    !is.na(region),
    # remove if in carehome on or before elig_date + 42 days (may need to reconsider carehome definition)
    is.na(longres_0_date),
    # remove if initialed end of life care on or before elig_date + 42 days
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

cat("#### clean vaccine data ####\n")
data_vaccine <- data_eligible %>%
  # calculate age based on JCVI group definition
  mutate(age = if_else(jcvi_group %in% c("10","11","12"), age_2, age_1)) %>%
  select(patient_id, jcvi_group, elig_date, age, region, starts_with("covid_vax")) %>%
  # rearrange to long format
  pivot_longer(cols = starts_with("covid")) %>%
  # extract brand name
  mutate(across(name, ~str_remove_all(., "covid_vax_|_date"))) %>%
  # remove if vax date missing
  filter(!is.na(value)) %>%
  # extract sequence of brand
  mutate(across(name, ~str_remove(.x, "_\\d"))) %>%
  # remove duplicates
  distinct() %>%
  # create dummy variable vax
  mutate(vax = TRUE) %>%
  # order each patient by value (vax date)
  arrange(patient_id, value) %>%
  # rearrange to wide with indicator variable for each brand (each row is a vax date)
  pivot_wider(names_from = name, values_from = vax) %>%
  mutate(across(c(moderna, disease, pfizer, az),
                ~if_else(is.na(.x), FALSE, .x))) %>%
  # only assign brand if that was the only brand recorded on that date (OK if disease code used on same date)
  mutate(brand = case_when(
    moderna & !(pfizer | az) ~ "moderna",
    pfizer & !(moderna | az) ~ "pfizer",
    az & !(moderna | pfizer) ~ "az",
    TRUE ~ "unknown"
  )) %>%
  # remove indicator variables
  select(-moderna, -disease, -pfizer, -az) %>%
  # rank the dates within each patient (shouldn't be any ties)
  group_by(patient_id) %>%
  mutate(dose = rank(value, ties.method = "min")) %>%
  ungroup() %>%
  # only keep 1st and 2nd doses
  filter(dose %in% 1:2) %>%
  # count number of each brand per patient
  group_by(patient_id, brand) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  # only keep if 2 doses, and 1st and 2nd dose either both pfizer or both az
  filter(n==2, brand %in% c("pfizer", "az")) %>%
  select(-n) %>%
  # rearrange to wide with a variable for each of 1st and 2nd dose
  pivot_wider(names_from = dose, values_from = value, names_prefix = "dose_") %>%
  # only keep if 2nd dose received [6,14) weeks after 1st dose
  filter(dose_1 + weeks(6) <= dose_2 & dose_2 < dose_1 + weeks(14))

cat("#### read elig_dates ####\n")
elig_dates <- readr::read_csv(here::here("output", "lib", "elig_dates.csv"))

# extract age range for each elig_date
elig_dates_age_range <- elig_dates %>%
  mutate(lower = str_extract(description, "age_. >= \\d{2}"),
         upper = str_extract(description, "age_. < \\d{2}")) %>%
  mutate(age_range = str_c(
    str_extract(lower, "\\d{2}"),
    " - ",
    str_extract(upper, "\\d{2}")
  )) %>%
  select(date, age_range)

source(here::here("analysis", "lib", "redaction_functions.R"))

cat("#### generate plot data ####\n")
# data to generate plots of second vax dates for all eligibility dates, stratified by region
out <- list()
i <- 1
for (plot_date in as.character(sort(unique(data_vaccine$elig_date)))) {
  
  ####
  
  data <- data_vaccine %>%
    filter(elig_date %in% as.Date(plot_date))
  
  # JCVI groups with the given plot_date
  jcvi_group_info <- data %>% 
    distinct(jcvi_group) %>%
    unlist() %>% unname() %>% sort() %>%
    str_c(., collapse = ", ")
  
  # age range for the given plot_date
  age_group_info <- elig_dates_age_range %>% 
    filter(date %in% as.Date(plot_date)) %>%
    select(age_range) %>%
    unlist() %>% unname()
  
  # sequence of dates for plot
  dates_seq <- seq(as.Date(plot_date) + weeks(6), 
                   as.Date(plot_date) + weeks(14) - days(1), 
                   1)
  
  # ensure full sequence of dates for each region:brand combo
  expanded_data <- tibble(
    region = character(),
    brand = character(),
    dose_2 = Date()
  )
  for (r in unique(data$region)) {
    for (v in unique(data$brand)) {
      expanded_data <- expanded_data %>%
        bind_rows(tibble(
          region = rep(r, each = length(dates_seq)),
          brand = rep(v, each = length(dates_seq)),
          dose_2 = dates_seq
        )) 
    }
  }
  
  # number of patients with 2nd dose on each date
  count_data <- data %>%
    group_by(region, brand, dose_2) %>%
    count() %>%
    ungroup() 
  
  # join expanded and count data
  plot_data <- expanded_data %>%
    left_join(count_data, by = c("region", "brand", "dose_2")) 
  
  plot_data_split <- plot_data %>%
    group_split(region, brand)
  
  # redact data within region:brand groups for plot
  out[[i]] <- bind_rows(
    lapply(
      plot_data_split, 
      function(x)
        x %>%
        # replace NAs with 0s
        mutate(across(n, ~if_else(is.na(.x), 0L, .x))) %>%
        # redact values < 5 (and next smallest if sum of redacted values <5)
        mutate(across(n, ~redactor2(.x))) %>%
        # replace redacted values with 5
        mutate(across(n, ~if_else(is.na(.x), 5L, .x)))
    )
  ) %>%
    mutate(jcvi_group = jcvi_group_info,
           age_group = age_group_info,
           elig_date = as.Date(plot_date, format = "%Y-%m-%d"))
  
  ####
  
  i <- i+1
}

# bind rows
plot_data_redacted <- bind_rows(out)

# save data for plotting
readr::write_csv(
  plot_data_redacted,
  file.path(data_dir,"plot_data_redacted.csv")
)
