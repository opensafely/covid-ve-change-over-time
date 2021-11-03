######################################

# This script plots the distribution of second vaccination dates within
# eligibility_date*regio strata

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
    
    elig_date = col_date(format="%Y-%m-%d"),
    jcvi_group = col_character(),
    
    # UK region
    region = col_character(),
    
    age_1 = col_integer(),
    age_2 = col_integer(),
    
    ## vaccination variables
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


# check format of elig_date
elig_date_test <- data_extract %>%
  select(elig_date) %>%
  filter(!is.na(elig_date) &
           str_detect(as.character(elig_date), "\\d{4}-\\d{2}-\\d{2}"))

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

cat("#### initial clean ####\n")
data_vaccines <- data_extract %>%
  filter(!(jcvi_group %in% "99"),
         !(elig_date %in% "2100-12-31")) %>%
  mutate(age = if_else(jcvi_group %in% c("10","11","12"), age_2, age_1)) %>%
  select(patient_id, jcvi_group, elig_date, age, region, starts_with("covid_vax")) %>%
  rename_with(.fn = ~str_remove_all(., "covid_vax_|_date"), .cols = starts_with("covid_vax")) %>%
  pivot_longer(
    cols = contains(c("disease", "pfizer", "az", "moderna")),
    names_to = c(".value", "n"),
    names_pattern = "(.+)_(.)"
  ) 

cat("#### filter vaccines to either pfizer or az ####\n")
data_keep <- data_vaccines %>%
  group_by(patient_id) %>%
  summarise(across(c("disease", "pfizer", "az", "moderna"), ~sum(!is.na(.)))) %>%
  ungroup() %>%
  # check with someone whether 'disease' should be used to check for errors
  select(-disease) %>%
  # remove all patients with a moderna vaccine
  filter(moderna == 0) %>%
  select(-moderna) %>%
  # remove if pfizer and az recorded for a single patient
  filter(!(pfizer > 0 & az > 0)) %>%
  select(patient_id)

cat("#### only plot if 2nd dose received [6,14) weeks after 1st dose ####\n")
data_remove <- data_extract %>%
  filter(
    !((covid_vax_pfizer_1_date + weeks(6) <= covid_vax_pfizer_2_date & 
       covid_vax_pfizer_2_date < covid_vax_pfizer_1_date + weeks(14)) |
      (covid_vax_az_1_date + weeks(6) <= covid_vax_az_2_date & 
         covid_vax_az_2_date < covid_vax_az_1_date + weeks(14)))
    ) %>%
  select(patient_id)

cat("#### identify 2nd dose dates ####\n")
data_vaccine_2 <- data_keep %>%
  left_join(data_vaccines, by = "patient_id") %>%
  anti_join(data_remove, by = "patient_id") %>%
  # only keep date of second vaccine
  filter(n %in% "2") %>%
  select(patient_id, elig_date, jcvi_group, age, region, pfizer, az) %>%
  pivot_longer(cols = c("pfizer", "az")) %>%
  filter(!is.na(value))


second_vax_dates_plot <- 
  function(
    plot_date = "2020-12-08"
  ) {
    
    jcvi_group_info <- data_vaccine_2 %>% 
      filter(elig_date %in% as.Date(plot_date)) %>%
      distinct(jcvi_group) %>%
      unlist() %>% unname() %>% sort() %>%
      str_c(., collapse = ", ")
    
    age_group_info <- data_vaccine_2 %>% 
      filter(elig_date %in% as.Date(plot_date)) %>%
      summarise(min = min(age), max = max(age)) %>%
      # transmute(range = str_c(min, " - ", max, " years")) %>%
      unlist() %>% unname()
    
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
    
    
    dates_seq <- seq(as.Date(plot_date) + weeks(6), 
                     as.Date(plot_date) + weeks(14) - days(1), 
                     1)
    
    expanded_data <- tibble(
      region = character(),
      name = character(),
      value = Date()
    )
    for (r in unique(data_vaccine_2$region)) {
      for (v in unique(data_vaccine_2$name)) {
        expanded_data <- expanded_data %>%
          bind_rows(tibble(
            region = rep(r, each = length(dates_seq)),
            name = rep(v, each = length(dates_seq)),
            value = dates_seq
            )) 
      }
    }
    
    count_data <- data_vaccine_2 %>%
      group_by(region, name, value) %>%
      count() %>%
      ungroup() 
    
    plot_data <- expanded_data %>%
      left_join(count_data, by = c("region", "name", "value")) %>%
      mutate(across(n, 
                    ~case_when(
                      .x < 100 ~ 100L,
                      is.na(.x) ~ 0L, 
                      TRUE ~ .x))) 
    
    x_breaks <- seq(as.Date(plot_date) + weeks(6),
                    as.Date(plot_date) + weeks(14),
                    14)
    
    
    plot_data %>%
      ggplot(aes(x = value, y = n, colour = name)) +
      geom_line() +
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
    
    ggsave(filename = here::here("output", "images", glue("second_vax_dates_{plot_date}.png")),
           width=18, height=14, units="cm")
    
  }

# plots of second vax dates for all eligibility dates, stratified by region
for (d in as.character(sort(unique(data_vaccine_2$elig_date)))) {
  second_vax_dates_plot(d)
}
