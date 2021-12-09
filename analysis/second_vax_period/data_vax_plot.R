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
  here::here("output", "vax", "data", "data_eligible_b.rds")
  )

data_vax_wide <- readr::read_rds(
  here::here("output", "vax", "data", "data_wide_vax_dates.rds")
  )

# second dose and brand for eligible individuals
data_2nd_dose <- data_eligible_b %>%
  left_join(data_vax_wide, by = "patient_id") %>%
  select(patient_id, elig_date, region_0, 
         dose_2 = covid_vax_2_date, brand = covid_vax_2_brand)
  
# function for creating plot data
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

# create list of data for each elig_date
data_vax_plot_list <- lapply(
  as.character(sort(unique(data_2nd_dose$elig_date))),
  function(x)
    try(generate_plot_data(plot_date = x))
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

thresholds <- data_vax_plot %>%
  group_by(elig_date, region_0, brand) %>%
  summarise(n = sum(n), .groups = "keep") %>%
  ungroup() %>%
  # threshold is the total number of individuals vaccinated in [+6wks, +16wks) divided by the number of days in the period (10*7)
  # i.e. start date is the first date that the 7-day average goes above 70-day average
  mutate(threshold = n/(10*7)) %>%
  select(-n)

# number of days in moving average
l <- 7 

# calculate moving averages for defining second vaccination periods and plotting
data_ma <- data_vax_plot %>%
  # calculate moving average number of individuals vaccinated for each elig_date:region_0:brand
  group_by(elig_date, region_0, brand) %>%
  mutate(moving_average = stats::filter(
    x = n,
    filter = rep(1/l, l),
    method = "convolution",
    sides = 2)) %>% # centred at day 4
  ungroup() %>%
  left_join(thresholds, by = c("elig_date", "region_0", "brand")) 

readr::write_rds(data_ma,
                 here::here("output", "second_vax_period", "data", "data_ma.rds"),
                 compress = "gz")

# second vaccination periods
second_vax_period_dates <- data_ma %>%
  # dates on which number of individuals vaccinated above threshold
  mutate(
    dose_2_above_threshold = if_else(
      moving_average > threshold,
      dose_2,
      NA_Date_
    )) %>%
  group_by(elig_date, region_0, brand) %>%
  mutate(
    # start the first date above threshold
    start_of_period = min(dose_2_above_threshold, na.rm = TRUE),
    # end the last date above threshold
    end_of_period = max(dose_2_above_threshold, na.rm = TRUE)
  ) %>%
  # only keep dates in period
  filter(dose_2 >= start_of_period, dose_2 <= end_of_period) %>%
  # count number of individuals vaccinated during period
  mutate(n_in_period = sum(n)) %>%
  ungroup() %>%
  # round to the closest 10 (so no need to redact, and reduce risk of secondary disclosure)
  mutate(across(n_in_period, ~ round(.x, -1))) %>%
  distinct(elig_date, region_0, brand, start_of_period, end_of_period, n_in_period)
# save a version to review and release
readr::write_csv(second_vax_period_dates,
                 here::here("output", "lib", "second_vax_period_dates.csv"))

# comparison dates for passing to study_definition_covs
comparison_dates <- second_vax_period_dates %>%
  # only keep if more than n_threshold individuals vaccinated with that brand in the second vaccination period
  filter(n_in_period > study_parameters$n_threshold) %>%
  # min start date / max end date for each elig_date/region, because cannot condition on vaccine brand in study_definition_covs
  group_by(elig_date, region_0, brand) %>%
  summarise(start_1_date = min(start_of_period) + days(14), 
            end_1_date = max(end_of_period) + days(14), 
            .groups = "keep") %>%
  ungroup() %>%
  mutate(condition = as.character(glue("(elig_date = {elig_date} AND region_0 = '{region_0}')"))) %>%
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
