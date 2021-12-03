######################################

# This script:
# - reads data_vax_plot.rds
# - identifies the second vaccination period
# - saves second_vax_period_dates.csv (the elig_date:region_0:brand specific dates)
# - saves start_dates.csv and end_dates.csv (the elig_date:region_0 specific dates to pass to study_definition_covs.py)
# - plots and saves the distribution of second vaccination dates

######################################

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  update_plots <- TRUE
  
} else {
  update_plots <- as.logical(args[[1]])
}

## setup
library(tidyverse)
library(lubridate)
library(glue)

# create folder for plots
images_dir <- here::here("output", "second_vax_period", "images")
dir.create(images_dir, showWarnings = FALSE, recursive=TRUE)

study_parameters <- readr::read_rds(here::here("output", "lib", "study_parameters.rds"))

# read data for plotting
data_vax_plot <- readr::read_rds(
  here::here("output", "second_vax_period", "data", "data_vax_plot.rds")
)

data_vax_plot <- data_vax_plot %>%
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

# parameters for plots
l <- 7 # number of days in moving average
n_threshold <- 100 # must have this number of individuals in second vax period to include the comparison for that elig_date/region/brand

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
  filter(n_in_period > n_threshold) %>%
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


if (update_plots) {
  
  # elig_dates info for plot titles
  elig_dates <- readr::read_csv(here::here("output", "lib", "elig_dates.csv")) %>%
    mutate(lower = str_extract(description, "age_. >= \\d{2}"),
           upper = str_extract(description, "age_. < \\d{2}")) %>%
    mutate(age_range = str_c(
      str_extract(lower, "\\d{2}"),
      " - ",
      str_extract(upper, "\\d{2}")
    )) %>%
    select(date, jcvi_groups, age_range)
  
  plot_2nd_vax_dates_fun <- function(
    data, 
    data_titles = elig_dates,
    plot_threshold = 5) {
    
    elig_date <- unique(data$elig_date)
    
    if (length(elig_date) != 1) stop("data$elig_date must contain one unique value")
    
    # JCVI groups and age range for plot title
    jcvi_groups <- data_titles %>%
      filter(date %in% elig_date) %>%
      select(jcvi_groups) %>% unlist() %>% unname()
    age_range <- data_titles %>%
      filter(date %in% elig_date) %>%
      select(age_range) %>% unlist() %>% unname()
    
    # plot title
    title_string <- glue("Patients eligible on {elig_date}")
    if (is.na(age_range)) {
      subtitle_string <- glue("JCVI group(s): {jcvi_groups}; Age range: whole group(s)")
    } else {
      subtitle_string <- glue("JCVI group(s): {jcvi_groups}; Age range: {age_range} years.")
    }
    
    # define breaks for x axis
    x_breaks <- seq(elig_date + weeks(6),
                    elig_date + weeks(16),
                    14) #14 days
    
    # plot histograms by region
    plot_by_region <- ggplot(NULL, aes(x = dose_2)) +
      # overlapping histogram for each brand, binwdith = 1 day
      geom_bar(data = data %>% filter(brand == "ChAdOx"), 
               aes(y = n, fill = "ChAdOx"),
               stat = "identity",  alpha = 0.5, width = 1) +
      geom_bar(data = data %>% filter(brand == "BNT162b2"), 
               aes(y = n, fill = "BNT162b2"), 
               stat = "identity", alpha = 0.5, width = 1) +
      # line for 7-day moving average for each brand
      geom_line(data = data %>% filter(brand == "ChAdOx") %>% filter(!is.na(moving_average)), 
                aes(y = moving_average, colour = "ChAdOx")) +
      geom_line(data = data %>% filter(brand == "BNT162b2") %>% filter(!is.na(moving_average)),  
                aes(y = moving_average, colour = "BNT162b2")) +
      # horizontal lines show the threshold above which the moving average must be fo the second vaccination period
      geom_hline(data = data %>% filter(brand == "ChAdOx"), 
                 aes(yintercept = threshold, colour = "ChAdOx"), 
                 linetype = "dashed") +
      geom_hline(data = data %>% filter(brand == "BNT162b2"), 
                 aes(yintercept = threshold, colour = "BNT162b2"), 
                 linetype = "dashed") +
      # facet by region
      facet_wrap(~ region_0, scales = "free_y") +
      scale_x_continuous(breaks = x_breaks,
                         labels = sapply(x_breaks, function(x) str_c(day(x), " ", month(x, label=TRUE)))) +
      scale_y_continuous(expand = expansion(mult = c(0,.05))) +
      scale_fill_discrete(guide = "none") +
      scale_colour_discrete(name = "brand") +
      labs(x = "date of second vaccination", y = "number of patients",
           title = title_string, subtitle = subtitle_string) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.x = element_text(size = 6),
            plot.caption = element_text(size = 10),
            plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.2, "cm")) +
      coord_cartesian(ylim = c(plot_threshold, NA))
    # caption:
    # X-axes restricted to 6 to 16 weeks after eligibility date.
    # Bars show the number of individuals who received a second dose of the given brand of vaccine on the given date.
    # Lines show the 7-day moving average (centred on day 4) of the numbers represented by the bars.
    
    # save the plot
    ggsave(plot_by_region,
           filename = file.path(images_dir, glue("plot_by_region_{elig_date}.png")),
           width=20, height=14, units="cm")
    
    
  }
  
  # generate and save plots
  lapply(data_ma %>% group_split(elig_date),
         function(x)
           try(plot_2nd_vax_dates_fun(x)))
  
}
