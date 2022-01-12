######################################

# This script:
# - reads data_eligible_b.rds and data_vax_wide.rds
# - saves data for plotting the distribution of 2nd vax dates across dates
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
                        labels = c("ChAdOx", "BNT162b2")))) %>%
  rename(n_brand = n) %>%
  group_by(jcvi_group, elig_date, region_0, dose_2) %>%
  mutate(n = sum(n_brand)) %>%
  ungroup()

readr::write_rds(data_vax_plot,
                 here::here("output", "second_vax_period", "data", "data_vax_plot.rds"),
                 compress = "gz")

# second vaccination periods
# number of days in cumulative sum
l <- 28 
second_vax_period_dates <- data_vax_plot %>%
  distinct(jcvi_group, elig_date, region_0, dose_2, n) %>%
  # calculate moving 28-day total number of individuals vaccinated for each elig_date:region_0:brand
  group_by(jcvi_group, elig_date, region_0) %>%
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
  mutate(
    # time between start of first comparison and last date of available data
    days_of_data = as.integer(as.Date(study_parameters$end_date) - start_of_period) + 14,
    # set n_comparisons based on days of available data
    n_comparisons = pmin(ceiling(days_of_data/28), study_parameters$max_comparisons)
    ) %>%
  select(-days_of_data)

brand_counts <- second_vax_period_dates %>%
  left_join(data_vax_plot,
            by = c("jcvi_group", "elig_date", "region_0")) %>%
  filter(
    start_of_period <= dose_2,
    dose_2 <= end_of_period
    ) %>%
  group_by(jcvi_group, elig_date, region_0, brand) %>%
  summarise(n = sum(n_brand), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(
    names_from = brand, values_from = n, names_prefix = "n_"
  )
  
second_vax_period_dates <- second_vax_period_dates %>%
  left_join(brand_counts,
            by = c("jcvi_group", "elig_date", "region_0"))  %>%
  select(jcvi_group, elig_date, region_0, n_ChAdOx, n_BNT162b2, cumulative_sum,
         start_of_period, end_of_period, n_comparisons)

# save for plotting
readr::write_rds(
  second_vax_period_dates,
  here::here("output", "lib", "second_vax_period_dates.rds"),
  compress = "gz")
# save to review and release, with cumulative sum rounded to nearest 10
capture.output(
  second_vax_period_dates %>% 
    mutate(across(cumulative_sum, ~ round(.x, -1))) %>% 
    kableExtra::kable("pipe"),
  file = here::here("output", "tables", "second_vax_period_dates.txt"),
  append=FALSE
)
