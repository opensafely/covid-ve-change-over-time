######################################

# This script:
# - reads data_plot.rds
# - plots the distribution of second vaccination dates within eligibility_date*region strata
# - saves the plots

######################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

# create folder for plots
images_dir <- here::here("output", "eda_index_dates", "images")
dir.create(images_dir, showWarnings = FALSE, recursive=TRUE)

# read data for plotting
data_plot <- readr::read_rds(
  here::here("output", "eda_index_dates", "data", "data_vax_plot.rds")
)#%>%
# mutate(across(n, ~.x*100))

# read elig_dates
elig_dates <- readr::read_csv(here::here("output", "lib", "elig_dates.csv")) %>%
  mutate(lower = str_extract(description, "age_. >= \\d{2}"),
         upper = str_extract(description, "age_. < \\d{2}")) %>%
  mutate(age_range = str_c(
    str_extract(lower, "\\d{2}"),
    " - ",
    str_extract(upper, "\\d{2}")
  )) %>%
  select(date, jcvi_groups, age_range)

# parameters for plots
l <- 7 # number of days in moving average
index_threshold <- 100 # threshold for starting and ending study period
plot_threshold <- 5 # mask counts <= plot_threshold 

second_vax_period_dates <- list()

i <- 1

cat("#### generate plots ####\n")
# plots of second vax dates for all eligibility dates, stratified by region
for (plot_date in as.character(sort(unique(data_plot$elig_date)))) {
  
  data <- data_plot %>%
    filter(elig_date %in% as.Date(plot_date)) %>%
    # change labels for plots    
    mutate(across(brand, 
                  ~factor(brand, 
                          levels = c("az", "pfizer"),
                          labels = c("ChAdOx", "BNT162b2")))) 
  
  jcvi_groups <- elig_dates %>%
    filter(date %in% as.Date(plot_date)) %>%
    select(jcvi_groups) %>% unlist() %>% unname()
  age_range <- elig_dates %>%
    filter(date %in% as.Date(plot_date)) %>%
    select(age_range) %>% unlist() %>% unname()

  # plot title
  title_string <- glue("Patients eligible on {plot_date}")
  if (is.na(age_range)) {
    subtitle_string <- glue("JCVI group(s): {jcvi_groups}; Age range: whole group(s)")
  } else {
    subtitle_string <- glue("JCVI group(s): {jcvi_groups}; Age range: {age_range} years.")
  }

  # define breaks for x axis
  x_breaks <- seq(as.Date(plot_date) + weeks(6),
                  as.Date(plot_date) + weeks(16),
                  14) #14 days
  
  data_ma <- data %>%
    group_by(region_0, brand) %>%
    mutate(moving_average = stats::filter(
      x = n,
      filter = rep(1/l, l),
      method = "convolution",
      sides = 2)) %>% 
    ungroup()
  

  
  second_vax_period_dates[[i]] <- data_ma %>%
    mutate(
      dose_2_above_threshold = if_else(
        moving_average > index_threshold,
        dose_2,
        NA_Date_
      )) %>%
    group_by(region_0, brand) %>%
    mutate(
      start_of_period = min(dose_2_above_threshold, na.rm = TRUE),
      end_of_period = max(dose_2_above_threshold, na.rm = TRUE)
      ) %>%
    filter(dose_2 >= start_of_period, dose_2 <= end_of_period) %>%
    mutate(n_in_period = sum(n)) %>%
    ungroup() %>%
    distinct(region_0, brand, start_of_period, end_of_period, n_in_period) %>%
    mutate(elig_date = as.Date(plot_date, format = "%Y-%m-%d"))
    
  plot_by_region <- ggplot(NULL, aes(x = dose_2)) +
    geom_bar(data = data_ma %>% filter(brand == "ChAdOx"), 
             aes(y = n, fill = "ChAdOx"),
             stat = "identity",  alpha = 0.5) +
    geom_bar(data = data_ma %>% filter(brand == "BNT162b2"), 
             aes(y = n, fill = "BNT162b2"), 
             stat = "identity", alpha = 0.5) +
    geom_line(data = data_ma %>% filter(brand == "ChAdOx"), 
              aes(y = moving_average, colour = "ChAdOx")) +
    geom_line(data = data_ma %>% filter(brand == "BNT162b2"), 
              aes(y = moving_average, colour = "BNT162b2")) +
    geom_hline(yintercept = index_threshold, 
               linetype = "dashed", colour = "grey") +
    facet_wrap(~ region_0, scales = "free_y") +
    scale_x_continuous(breaks = x_breaks,
                       labels = sapply(x_breaks, function(x) str_c(day(x), " ", month(x, label=TRUE)))) +
    scale_y_continuous(expand = expansion(mult = c(0,.05))) +
    scale_fill_discrete(guide = FALSE) +
    scale_colour_discrete(name = "brand") +
    labs(x = "date of second vaccination", y = "number of patients",
         title = title_string, subtitle = subtitle_string,
         caption = caption_string) +
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
         filename = file.path(images_dir, glue("plot_by_region_{plot_date}.png")),
         width=20, height=14, units="cm")
    
  i <- i+1
}

# save the start and end dates of the second vax period 
second_vax_period_dates <- bind_rows(second_vax_period_dates)
readr::write_csv(second_vax_period_dates,
                 here::here("output", "lib", "second_vax_period_dates.csv"))
