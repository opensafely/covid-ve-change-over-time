######################################

# This script:
# - reads data_2nd_dose.rds
# - generates and saves data_plot for plotting the distribution of 2nd vax dates
######################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

data_eligible_b <- readr::read_rds(
  here::here("output", "vax", "data", "data_eligible_b.rds")
  )

data_vax_wide <- readr::read_rds(
  here::here("output", "vax", "data", "data_wide_vax_dates.rds")
  )

data_2nd_dose <- data_eligible_b %>%
  left_join(data_vax_wide, by = "patient_id") %>%
  select(patient_id, elig_date, region_0, 
         dose_2 = covid_vax_2_date, brand = covid_vax_2_brand)
  

generate_plot_data <- function(data = data_2nd_dose, plot_date) {
  
  # check date in correct format
  d <- try(as.Date(plot_date, format="%Y-%m-%d"))
  if("try-error" %in% class(d) || is.na(d)) stop("plot_date is not in %Y-%m-%d format.")
  
  data <- data %>%
    filter(elig_date %in% as.Date(plot_date))
  
  if (nrow(data)==0) stop("No samples for the given plot_date.")
  
  # sequence of dates for plot
  dates_seq <- seq(as.Date(plot_date) + weeks(6), 
                   as.Date(plot_date) + weeks(16) - days(1), 
                   1)
  
  # ensure full sequence of dates for each region:brand combo
  # so that the moving averages calculated in 'plot_2nd_vax_dates.R' include the zero counts
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
  out <- expanded_data %>%
    left_join(count_data, by = c("region_0", "brand", "dose_2")) %>%
    mutate(across(n, ~if_else(is.na(.x), 0L, .x))) %>%
    mutate(elig_date = as.Date(plot_date, format = "%Y-%m-%d"))
  
  return(out)
  
}

# create list of plot_data for each elig_date
data_vax_plot_list <- lapply(
  as.character(sort(unique(data_2nd_dose$elig_date))),
  function(x)
    try(generate_plot_data(plot_date = x))
)

data_vax_plot <- bind_rows(data_vax_plot_list[sapply(data_vax_plot_list, is_tibble)])

# save data for plotting
readr::write_rds(
  data_vax_plot,
  here::here("output", "second_vax_period", "data", "data_vax_plot.rds"),
  compress = "gz"
)
