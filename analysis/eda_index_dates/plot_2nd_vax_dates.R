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
)

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
  
  caption <- "X-axes show 6 to 16 weeks after eligibility date. Bars are stacked by vaccine brand."

  # plot the data
  plot_by_region <- data %>%
    # calculate approx. total within region:brand group (approx as post-redaction)
    group_by(region_0) %>%
    mutate(n_region = scales::comma(sum(n), accuracy = 1)) %>%
    ungroup() %>%
    mutate(across(region_0,
                  ~str_replace(.x,
                               "Yorkshire and The Humber",
                               "Yorkshire & Humber"))) %>%
    mutate(across(region_0, ~glue("{region_0} (n={n_region})"))) %>%
    ggplot(aes(x = dose_2, y = n, fill = brand)) +
    geom_bar(stat = "identity", position = "stack", width=1, alpha=0.9) +
    facet_wrap(~ region_0, scales = "free_y") +
    labs(x = "date of second vaccination", y = "number of patients",
         title = title_string, subtitle = subtitle_string,
         caption = str_wrap(caption, 120)) +
    scale_fill_discrete(name = "vaccine brand") +
    scale_x_continuous(breaks = x_breaks,
                       labels = sapply(x_breaks, function(x) str_c(day(x), " ", month(x, label=TRUE)))) +
    scale_y_continuous(expand = expansion(mult = c(0,.05))) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 6),
          plot.caption = element_text(size = 10),
          plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.2, "cm")) +
    coord_cartesian(ylim = c(plot_threshold, NA))

  # save the plot
  ggsave(plot_by_region,
         filename = file.path(images_dir, glue("plot_by_region_{plot_date}.png")),
         width=20, height=14, units="cm")
  
  data_plot_all <- data %>%
    group_by(brand, dose_2) %>%
    summarise(n_brand = sum(n), .groups = "keep") %>%
    ungroup() %>%
    group_by(brand) %>%
    mutate(moving_average = stats::filter(
      x = n_brand,
      filter = rep(1/l, l),
      method = "convolution",
      sides = 2)) 
  
  second_vax_period_dates[[i]] <- data_plot_all %>%
    # distinct(brand, dose_2, moving_average) %>%
    filter(moving_average>index_threshold) %>%
    group_by(brand) %>%
    summarise(across(dose_2, list(min=min,max=max)), .groups = "keep") %>%
    ungroup() %>%
    mutate(elig_date = as.Date(plot_date, format = "%Y-%m-%d"))
  
  plot_all <- data_plot_all %>%
    left_join(second_vax_period_dates[[i]], 
              by = "brand") %>%
    mutate(above_threshold = if_else(
      moving_average > index_threshold,
      dose_2,
      NA_Date_)) %>%
    group_by(brand) %>%
    mutate(start = min(above_threshold, na.rm=TRUE),
           end = max(above_threshold, na.rm=TRUE)) %>%
    mutate(fill = if_else(dose_2 >= start & dose_2 <= end,
                          TRUE, FALSE)) %>%
    # mutate(keep_brand = any(fill)) %>%
    ungroup() %>%
    # filter(keep_brand) %>%
    droplevels() %>%
    ggplot() +
    geom_bar(aes(x = dose_2, y = n_brand, fill = fill),
             stat = "identity", position = "stack", width=1, alpha=0.9) +
    geom_line(aes(x = dose_2, y = moving_average)) +
    geom_hline(yintercept = index_threshold, linetype = "dashed") +
    facet_wrap(~brand, scales = "free_y", nrow=2) +
    labs(x = "date of second vaccination", y = "number of patients",
         title = title_string, subtitle = subtitle_string,
         caption = str_wrap(glue("Solid lines show {l}-day moving average of number of patients. The start of the \"second vaccination period\" is the date at which the solid line first crosses the dashed line (n = {index_threshold}), and the end of the study period is the last point at which the solid line crosses the dashed line. X-axes show 6 to 16 weeks after eligibility date."), 120)) +
    scale_fill_manual(breaks = c(FALSE, TRUE), 
                      values=c("#bdbdbd", "#31a354"),
                      name = "Included in \"second vaccination period\":") +
    scale_x_continuous(breaks = x_breaks,
                       labels = sapply(x_breaks, function(x) str_c(day(x), " ", month(x, label=TRUE)))) +
    scale_y_continuous(expand = expansion(mult = c(0,.05))) +
    guides(fill = guide_legend(title.position = "top")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 6),
          plot.caption = element_text(size = 10),
          plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.2, "cm")) +
    coord_cartesian(ylim = c(plot_threshold, NA))
  
  # save the plot
  ggsave(plot_all,
         filename = file.path(images_dir, glue("plot_all_{plot_date}.png")),
         width=20, height=14, units="cm")
    
  i <- i+1
}

# save the start and end dates of the second vax period 
second_vax_period_dates <- bind_rows(second_vax_period_dates)
readr::write_csv(second_vax_period_dates,
                 here::here("output", "lib", "second_vax_period_dates.csv"))
