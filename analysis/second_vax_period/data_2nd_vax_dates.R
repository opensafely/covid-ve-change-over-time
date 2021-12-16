######################################

# This script:
# - reads data_eligible_b.rds and data_vax_wide.rds
# - saves data_ma.rds for plotting the distribution of 2nd vax dates across dates
# - identifies the second vaccination period
# - saves second_vax_period_dates.csv (the elig_date:region_0:brand specific dates)
# - saves start_dates.csv and end_dates.csv (the elig_date:region_0 specific dates to pass to study_definition_covs.py)

######################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

# create folder for data
images_dir <- here::here("output", "second_vax_period", "data")
dir.create(images_dir, showWarnings = FALSE, recursive=TRUE)

study_parameters <- readr::read_rds(here::here("output", "lib", "study_parameters.rds"))

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_b <- readr::read_rds(
  here::here("output", "data", "data_eligible_b.rds")
  )

data_vax_wide <- readr::read_rds(
  here::here("output", "data", "data_wide_vax_dates.rds")
  )

# second dose and brand for eligible individuals
data_2nd_dose <- data_eligible_b %>%
  left_join(data_vax_wide, by = "patient_id") %>%
  select(patient_id, jcvi_group, elig_date, region_0, 
         dose_2 = covid_vax_2_date, brand = covid_vax_2_brand) %>%
  group_split(jcvi_group, elig_date)

  
# function for creating plot data
generate_plot_data <- function(.data) {
  
  group <- unique(.data$jcvi_group)
  plot_date <- unique(.data$elig_date)
  
  if (nrow(.data)==0) stop(".data is an empty tibble.")
  
  # sequence of dates for plot
  dates_seq <- seq(as.Date(plot_date) + weeks(6), 
                   as.Date(plot_date) + weeks(20) - days(1), 
                   1)
  
  # ensure full sequence of dates for each region:brand combo
  # so that the moving averages calculated in 'plot_2nd_vax_dates.R' include the zero counts
  expanded_data <- tibble(
    region_0 = character(),
    brand = character(),
    dose_2 = Date()
  )
  for (r in unique(.data$region_0)) {
    for (v in unique(.data$brand)) {
      expanded_data <- expanded_data %>%
        bind_rows(tibble(
          region_0 = rep(r, each = length(dates_seq)),
          brand = rep(v, each = length(dates_seq)),
          dose_2 = dates_seq
        ))
    }
  }
  
  # number of patients with 2nd dose on each date
  count_data <- .data %>%
    group_by(region_0, brand, dose_2) %>%
    count() %>%
    ungroup() 
  
  # join expanded and count data
  out <- expanded_data %>%
    left_join(count_data, by = c("region_0", "brand", "dose_2")) %>%
    mutate(across(n, ~if_else(is.na(.x), 0L, .x))) %>%
    mutate(
      jcvi_group = group,
      elig_date = as.Date(plot_date, format = "%Y-%m-%d"))
  
  return(out)
  
}

# create list of data for each elig_date
data_vax_plot_list <- lapply(
  data_2nd_dose,
  function(x)
    try(generate_plot_data(x))
)

# bind list into one tibble
data_vax_plot <- bind_rows(
  data_vax_plot_list[sapply(data_vax_plot_list, is_tibble)]
  ) %>%
  # change labels for plots    
  mutate(across(brand, 
                ~factor(brand, 
                        levels = c("az", "pfizer"),
                        labels = c("ChAdOx", "BNT162b2")))) 

readr::write_rds(data_vax_plot,
                 here::here("output", "second_vax_period", "data", "data_vax_plot.rds"),
                 compress = "gz")

# second vaccination periods
# number of days in cumulative sum
l <- 28 
second_vax_period_dates <- data_vax_plot %>%
  # calculate moving average number of individuals vaccinated for each elig_date:region_0:brand
  group_by(jcvi_group, elig_date, region_0, brand) %>%
  arrange(dose_2, .by_group = TRUE) %>%
  mutate(
    
    cumulative_sum = stats::filter(
      x = n,
      filter = rep(1, l),
      method = "convolution",
      sides = 1),
    
    end_of_period = if_else(
      cumulative_sum == max(cumulative_sum, na.rm = TRUE),
      dose_2,
      as.Date(NA_character_)),
    
    start_of_period = end_of_period - days(27)
    
  ) %>% 
  # only keep rows where cumulative_sum == max(cumulative_sum, na.rm = TRUE)
  filter(!is.na(end_of_period)) %>%
  # in case there are multiple dates with max(cumulative_sum),
  # take the first date with max(cumulative_sum)
  summarise(across(c(cumulative_sum, end_of_period, start_of_period), 
                   min, na.rm = TRUE),
            .groups = "keep") %>%
  ungroup() %>%
  distinct(jcvi_group, brand, elig_date, region_0, 
           start_of_period, end_of_period, cumulative_sum)

# save for plotting
readr::write_rds(
  second_vax_period_dates,
  here::here("output", "lib", "second_vax_period_dates.rds"),
  compress = "gz")
# save to review and release, with cumulative sum rounded to nearest 10
readr::write_csv(
  second_vax_period_dates %>% mutate(across(cumulative_sum, ~ round(.x, -1))),
  here::here("output", "lib", "second_vax_period_dates.csv"))

# comparison dates for passing to study_definition_covs
comparison_dates <- second_vax_period_dates %>%
  # only keep if more than n_threshold individuals vaccinated in the jcvi_group:elig_date:region:brand period
  filter(cumulative_sum > study_parameters$n_threshold) %>%
  # min start date / max end date for each elig_date/region, because cannot condition on vaccine brand in study_definition_covs
  group_by(jcvi_group, elig_date, region_0, brand) %>%
  summarise(start_1_date = min(start_of_period) + days(14), 
            end_1_date = max(end_of_period) + days(14), 
            .groups = "keep") %>%
  ungroup() %>%
  mutate(condition = as.character(glue("(jcvi_group = '{jcvi_group}' AND elig_date = {elig_date} AND region_0 = '{region_0}')"))) %>%
  select(start_1_date, end_1_date, condition) %>%
  add_row(start_1_date = as.Date("2100-01-01"), 
          end_1_date = as.Date("2100-12-31"), 
          condition = "DEFAULT")

start_dates <- comparison_dates %>%
  select(-end_1_date) %>%
  arrange(start_1_date) %>%
  group_by(start_1_date) %>%
  summarise(condition = str_c(condition, collapse  = " OR "), .groups = "keep") %>%
  ungroup()

end_dates <- comparison_dates %>%
  select(-start_1_date) %>%
  arrange(end_1_date) %>%
  group_by(end_1_date) %>%
  summarise(condition = str_c(condition, collapse  = " OR "), .groups = "keep") %>%
  ungroup()

# save for passing to study_definition_covs.py
readr::write_csv(start_dates,
                 here::here("output", "lib", "start_dates.csv"))

readr::write_csv(end_dates,
                 here::here("output", "lib", "end_dates.csv"))
