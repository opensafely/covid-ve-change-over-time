######################################

# This script:
# - reads data_2nd_dose.rds
# - generates and saves data_plot for plotting the distribution of 2nd vax dates
######################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

data_2nd_dose <- readr::read_rds(
  here::here("output", "eda_index_dates", "data", "data_2nd_dose.rds"))

# data to generate plots of second vax dates for all eligibility dates, stratified by region
out <- list()
i <- 1
for (plot_date in as.character(sort(unique(data_2nd_dose$elig_date)))) {
  
  data <- data_2nd_dose %>%
    filter(elig_date %in% as.Date(plot_date))
  
  # sequence of dates for plot
  dates_seq <- seq(as.Date(plot_date) + weeks(6), 
                   as.Date(plot_date) + weeks(16) - days(1), 
                   1)
  
  # ensure full sequence of dates for each region:brand combo
  expanded_data <- tibble(
    region_0 = character(),
    brand = character(),
    dose_2 = Date()
  )
  for (r in unique(data$region_0)) {
    for (v in unique(data$brand)) {
      expanded_data <- expanded_data %>%
        bind_rows(tibble(
          region_0 = rep(r, each = length(dates_seq)),
          brand = rep(v, each = length(dates_seq)),
          dose_2 = dates_seq
        )) 
    }
  }
  
  # number of patients with 2nd dose on each date
  count_data <- data %>%
    group_by(region_0, brand, dose_2) %>%
    count() %>%
    ungroup() 
  
  # join expanded and count data
  out[[i]] <- expanded_data %>%
    left_join(count_data, by = c("region_0", "brand", "dose_2")) %>%
    mutate(across(n, ~if_else(is.na(.x), 0L, .x))) %>%
    mutate(elig_date = as.Date(plot_date, format = "%Y-%m-%d"))
  
  i <- i+1
}

# bind rows
data_plot <- bind_rows(out)

# save data for plotting
readr::write_rds(
  data_plot,
  here::here("output", "eda_index_dates", "data", "data_plot.rds")
)