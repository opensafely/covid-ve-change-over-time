######################################

# This script:
# - reads the extracted data
# - processes the extracted data
# - cleans the vaccination data to identify second doses of pfizer and az
# - plots the distribution of second vaccination dates within eligibility_date*region strata

######################################

## setup
library(tidyverse)
library(lubridate)
library(readr)
library(glue)

## create folders for outputs
dir.create(here::here("output", "data"), showWarnings = FALSE, recursive=TRUE)
dir.create(here::here("output", "tables"), showWarnings = FALSE, recursive=TRUE)
dir.create(here::here("output", "images"), showWarnings = FALSE, recursive=TRUE)

## import dates
dates <- readr::read_rds(here::here("analysis", "lib", "study_dates.rds"))

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
    
    ## Eligibiloty dates
    elig_date = col_date(format="%Y-%m-%d"),
    jcvi_group = col_character(),
    
    ## UK region
    region = col_character(),
    
    ## Age
    age_1 = col_integer(),
    age_2 = col_integer(),
    
    ## Variables for applying exclusion criteria
    positive_test_0_date = col_date(format="%Y-%m-%d"),
    primary_care_covid_case_0_date = col_date(format="%Y-%m-%d"),
    primary_care_suspected_covid_0_date = col_date(format="%Y-%m-%d"),
    covidadmitted_0_date = col_date(format="%Y-%m-%d"),
    longres_0_date = col_date(format="%Y-%m-%d"),
    endoflife_0_date = col_date(format="%Y-%m-%d"),
    midazolam_0_date = col_date(format="%Y-%m-%d"),
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
  jcvi_groups <- readr::read_csv(here::here("analysis", "lib", "jcvi_groups.csv"))
  
  elig_dates <- readr::read_csv(here::here("analysis", "lib", "elig_dates.csv"))
  
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
  data_extract <- eval(parse(text = glue("data_extract %>% mutate(elig_date = as.Date(case_when({elig_date_cases})), format=\'%Y-%m-%d\')")))
  
  data_extract <- data_extract %>%
    mutate(r1 = round(rnorm(nrow(.), mean = 7, sd = 7)),
           r2 = round(rnorm(nrow(.), mean = 7, sd = 7)),
           r3 = round(rnorm(nrow(.), mean = 8*7, sd = 7)),
           r4 = round(rnorm(nrow(.), mean = 8*7, sd = 7))) %>%
    mutate(across(covid_vax_pfizer_1_date, ~if_else(!is.na(.x), elig_date + days(r1),.x))) %>%
    mutate(across(covid_vax_az_1_date, ~if_else(!is.na(.x), elig_date + days(r2),.x))) %>%
    mutate(covid_vax_pfizer_2_date = covid_vax_pfizer_1_date + days(r3),
           covid_vax_az_2_date = covid_vax_az_1_date + days(r4))

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
      TRUE ~ "Missing"
    ),
  ) %>%
  # apply exclusion criteria
  filter(
    # remove if any missing data for key variables
    !(ethnicity %in% "Missing"),
    !(sex %in% ""),
    !(imd %in% ""),
    !(region %in% ""),
    # remove if in carehome on or before elig_date + 42 days (may need to reconsider carehome definition)
    is.na(longres_0_date),
    # remove if initialed end of life care on or before elig_date + 42 days
    is.na(endoflife_0_date),
    is.na(midazolam_0_date),
    # remove if evidence of covid infection on or before elig_date + 42 days
    is.na(positive_test_0_date),
    is.na(primary_care_covid_case_0_date),
    is.na(primary_care_suspected_covid_0_date),
    is.na(covidadmitted_0_date)
  )
  

cat("#### clean vaccine data ####\n")
data_vaccine <- data_processed %>%
  # remove dummy groups and dates
  filter(!(jcvi_group %in% "99"),
         !(elig_date %in% "2100-12-31")) %>%
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

# function for plotting distribution of 2nd vax dates
second_vax_dates_plot <- 
  function(
    plot_date = "2020-12-08"
  ) {
    
    # JCVI groups with the given plot_date
    jcvi_group_info <- data_vaccine %>% 
      filter(elig_date %in% as.Date(plot_date)) %>%
      distinct(jcvi_group) %>%
      unlist() %>% unname() %>% sort() %>%
      str_c(., collapse = ", ")
    
    # age range for the given plot_date
    age_group_info <- data_vaccine %>% 
      filter(elig_date %in% as.Date(plot_date)) %>%
      summarise(min = min(age), max = max(age)) %>%
      # transmute(range = str_c(min, " - ", max, " years")) %>%
      unlist() %>% unname()
    
    # plot title
    title_string <- glue("Patients eligible on {plot_date}")
    if (age_group_info[1]==age_group_info[2]) {
      age_group_info <- glue("{age_group_info[1]} years")
      subtitle_string <- glue("JCVI group(s): {jcvi_group_info}; Age: {age_group_info}.")
    } else if (age_group_info[2] > 90) {
      age_group_info <- str_c(age_group_info[1], "+ years")
      subtitle_string <- glue("JCVI group(s): {jcvi_group_info}; Age range: {age_group_info}.")
    } else {
      age_group_info <- str_c(age_group_info[1], " - ", age_group_info[2], " years")
      subtitle_string <- glue("JCVI group(s): {jcvi_group_info}; Age range: {age_group_info}.")
    }
    
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
    for (r in unique(data_vaccine$region)) {
      for (v in unique(data_vaccine$brand)) {
        expanded_data <- expanded_data %>%
          bind_rows(tibble(
            region = rep(r, each = length(dates_seq)),
            brand = rep(v, each = length(dates_seq)),
            dose_2 = dates_seq
            )) 
      }
    }
    
    # number of patients with 2nd dose on each date
    count_data <- data_vaccine %>%
      group_by(region, brand, dose_2) %>%
      count() %>%
      ungroup() 
    
    # join and mask dates on which <10 patients received their 2nd dose
    plot_data <- expanded_data %>%
      left_join(count_data, by = c("region", "brand", "dose_2")) %>%
      mutate(across(n, 
                    ~case_when(
                      .x < 10 ~ 10L,
                      is.na(.x) ~ 0L, 
                      TRUE ~ .x))) 
    
    # define breaks for x axis
    x_breaks <- seq(as.Date(plot_date) + weeks(6),
                    as.Date(plot_date) + weeks(14),
                    14)
    
    # plot the data
    plot_data %>%
      ggplot(aes(x = dose_2, y = n, colour = brand)) +
      geom_line() +
      # line at elig_date + 10 weeks, as this is potentially going to be time_zero for comparisons
      geom_vline(xintercept = as.Date(plot_date, format = "%Y-%m-%d") + weeks(10),
                 linetype = "dashed") +
      facet_wrap(~ region, scales = "free_y") +
      labs(x = "date of second vaccination", y = "number of patients",
           title = title_string, subtitle = subtitle_string) +
      scale_color_discrete(name = "vaccine") +
      scale_x_continuous(breaks = x_breaks,
                         labels = sapply(x_breaks, function(x) str_c(day(x), " ", month(x, label=TRUE)))) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            axis.text.x = element_text(size = 6),
            plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.2, "cm")) +
      # to avoid displaying low numbers of patients
      coord_cartesian(ylim = c(10, NA))
    
    # save the plot
    ggsave(filename = here::here("output", "images", glue("second_vax_dates_{plot_date}.png")),
           width=18, height=14, units="cm")
    
  }

cat("#### generate plots ####\n")
# plots of second vax dates for all eligibility dates, stratified by region
for (d in as.character(sort(unique(data_vaccine$elig_date)))) {
  second_vax_dates_plot(d)
}
