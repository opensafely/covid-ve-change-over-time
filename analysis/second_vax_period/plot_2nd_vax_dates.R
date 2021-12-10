######################################

# This script:
# - plots and saves the distribution of second vaccination dates

######################################

library(tidyverse)
library(glue)
library(lubridate)

# create folder for plots
images_dir <- here::here("output", "second_vax_period", "images")
dir.create(images_dir, showWarnings = FALSE, recursive=TRUE)

data_ma <- readr::read_rds(here::here("output", "second_vax_period", "data", "data_ma.rds"))

# elig_dates info for plot titles
group_age_ranges <- readr::read_rds(
  here::here("output", "lib", "group_age_ranges.rds"))


plot_2nd_vax_dates_fun <- function(
  data, 
  subtitle_string = group_age_ranges,
  plot_threshold = 5) {
  
  jcvi_group <- unique(data$jcvi_group)
  elig_date <- unique(data$elig_date)
  
  # plot title
  title_string <- glue("JCVI group {jcvi_group}; eligible from {elig_date}")
  
  # age range for plot title
  subtitle_string <- str_c(
    "Age range: ",
    subtitle_string$age_range[
      subtitle_string$elig_date == elig_date &
        subtitle_string$jcvi_group == jcvi_group
    ],
    " years"
  )
  
  # define breaks for x axis
  x_breaks <- seq(elig_date + weeks(6),
                  elig_date + weeks(20),
                  28) #28 days
  
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
         filename = file.path(images_dir, glue("plot_by_region_{jcvi_group}_{elig_date}.png")),
         width=20, height=14, units="cm")
  
}

# generate and save plots
lapply(data_ma %>% group_split(jcvi_group, elig_date),
       function(x)
         try(plot_2nd_vax_dates_fun(data = x)))

