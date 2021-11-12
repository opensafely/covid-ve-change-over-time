######################################

# This script:
# - reads plot_data_extracted.csv
# - plots the distribution of second vaccination dates within eligibility_date*region strata
# - saves the plots

######################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

# create folder for plots
images_dir <- here::here("output", "explore_2nd_vax_dates", "images")
dir.create(images_dir, showWarnings = FALSE, recursive=TRUE)

# read data for plotting
plot_data_redacted <- readr::read_csv(here::here("output", "explore_2nd_vax_dates", "data", "plot_data_redacted.csv"))

cat("#### generate plots ####\n")
# plots of second vax dates for all eligibility dates, stratified by region
for (plot_date in as.character(sort(unique(plot_data_redacted$elig_date)))[1]) {
  
  data <- plot_data_redacted %>%
    filter(elig_date %in% as.Date(plot_date))
  
  jcvi_group_info <- unique(data$jcvi_group)
  age_group_info <- unique(data$age_group)
  
  # plot title
  title_string <- glue("Patients eligible on {plot_date}")
  if (is.na(age_group_info)) {
    subtitle_string <- glue("JCVI group(s): {jcvi_group_info}; Age range: whole group(s)")
  } else {
    subtitle_string <- glue("JCVI group(s): {jcvi_group_info}; Age range: {age_group_info} years.")
  }
  
  # define breaks for x axis
  x_breaks <- seq(as.Date(plot_date) + weeks(6),
                  as.Date(plot_date) + weeks(14),
                  14)
  
  # plot the data
  data %>%
    # calculate approx. total within region:brand group (approx as post-redaction)
    group_by(region) %>%
    mutate(n_region = scales::comma(sum(n), accuracy = 1)) %>%
    ungroup() %>%
    mutate(across(region,
                  ~str_replace(.x,
                               "Yorkshire and The Humber",
                               "Yorkshire & Humber"))) %>%
    mutate(across(region, ~glue("{region} (n={n_region})"))) %>%
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
          plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.2, "cm"))
  
  # save the plot
  ggsave(filename = file.path(images_dir, glue("second_vax_dates_{plot_date}.png")),
         width=20, height=14, units="cm")
}
