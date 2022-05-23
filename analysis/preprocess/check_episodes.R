################################################################################
# load libraries
library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(glue)

################################################################################
# read metadata
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

################################################################################
# load functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
# load data
data_episodes <- readr::read_rds(
  here::here("output", "data", "data_episodes.rds"))

################################################################################
# create output directory
fs::dir_create(here::here("output", "eda"))

################################################################################
# bar plot of episode triggers

# define time periods for bar plot
weeks <- seq(
  as.Date("2020-03-01"), 
  as.Date(study_parameters$end_date) + weeks(2), 
  by = 14
  )

# prepare data for bar plot
data_bar <- data_episodes %>%
  select(-episode, -episode_end_date) %>%
  mutate(across(
    c(starts_with("covid_primary_care"), covidadmitted_date, postest_date, coviddeath_date), 
    ~ .x == episode_start_date & !is.na(.x))) %>%
  pivot_longer(
    cols = c(starts_with("covid_primary_care"), covidadmitted_date, postest_date, coviddeath_date)
    ) %>%
  filter(value) %>%
  mutate(across(
    name,
    ~ factor(
      .x, 
      levels = str_c(c("postest", "covidadmitted", "covid_primary_care_positive_test", "covid_primary_care_code", "covid_primary_care_sequalae", "coviddeath"), "_date"),
      labels = c("Positive test (SGSS)", "COVID-19 hospitalisation", "Positive test (primary care)", "Diagnosis (primary care)", "Sequalae (primary care)", "COVID-19 death")
    ))) %>%
  # if there are multiple events on the episode start date, keep in the order of the levels of the name factor, defined above
  arrange(patient_id, episode_start_date, name) %>%
  distinct(patient_id, episode_start_date, .keep_all = TRUE) %>%
  filter(
    as.Date("2020-03-01") <= episode_start_date,
    episode_start_date <= as.Date(study_parameters$end_date)
    ) %>%
  mutate(across(episode_start_date, ~cut(.x, breaks = weeks, right = FALSE))) %>%
  droplevels() %>%
  group_by(episode_start_date, name) %>%
  count() %>%
  ungroup(name) %>%
  mutate(across(n, ~ceiling_any(.x, to=7))) %>%
  mutate(
    id = row_number(),
    total = sum(n),
    prop = n/total
    ) %>%
  ungroup() %>%
  mutate(across(total, 
                ~if_else(
                  id==1,
                  scales::comma(.x, accuracy = 1),
                  NA_character_
                  )))

# cutoff the plot to make the less common events more visable
# x_upper <- 1-min(data_bar$prop[data_bar$name == "Positive test (SGSS)"])
x_upper <- 0.4

max_width <- max(nchar(data_bar$total), na.rm = TRUE)

# create bar plot
data_bar %>%
  mutate(across(total, ~str_pad(.x, width = max_width, side = "left", pad = " "))) %>%
  ggplot(aes(y = reorder(episode_start_date, -as.integer(episode_start_date)))) +
  geom_bar(aes(x = prop, fill = name), stat = "identity", width = 1) +
  geom_text(aes(x = x_upper-0.02, label = total), size = 2.5) +
  scale_y_discrete(
    name = "Episode start date"
  ) +
  scale_x_continuous(
    name = NULL,
    labels = scales::percent_format()
    ) +
  coord_cartesian(xlim = c(0, x_upper)) +
  guides(
    fill = guide_legend(
      title = "Episode trigger",
      nrow=3
      )) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(margin = margin(c(t=10,r=0,b=0,l=0))),
    legend.position = "bottom"
  )

ggsave(
  filename = here::here("output", "eda", "episode_triggers.png"),
                        width = 18, height = 26, units = "cm")

################################################################################
# scatter plot of episode length vs episode start date 

episode_length_plot <- function(trigger = NULL) {
  
  if (!is.null(trigger)) {
    
    data_tile <- data_episodes %>%
      filter(episode_start_date == !! sym(trigger))
    
  } else {
    data_tile <- data_episodes
    trigger <- "any"
  }
  
  x_tilewidth <- 7
  y_tilewidth <- 10
  
  data_tile <- data_tile %>%
    mutate(episode_length = as.numeric(episode_end_date - episode_start_date)) %>%
    mutate(across(episode_start_date, ~ cut(.x, breaks = seq(min(.x), max(.x) + days(x_tilewidth), by = x_tilewidth), include.lowest = TRUE))) %>%
    mutate(across(episode_length, ~ cut(.x, breaks = seq(0, max(.x) + y_tilewidth, by = y_tilewidth), include.lowest = TRUE))) 
  
  x_breaks <- levels(data_tile$episode_start_date)
  x_breaks <- x_breaks[seq(10,length(x_breaks),10)]
  y_breaks <- levels(data_tile$episode_length)
  y_breaks <- y_breaks[seq(11,length(y_breaks),10)]
  y_breaks_labs <- str_extract(y_breaks, "\\d+")
  
  
  data_tile %>%
    group_by(episode_start_date, episode_length) %>%
    count() %>%
    ungroup() %>%
    mutate(across(n, ~ceiling_any(.x, to = 7))) %>%
    ggplot(aes(x = episode_start_date, y = episode_length, fill = n)) +
    geom_tile() +
    scale_x_discrete(
      name = "Episode start date",
      breaks = x_breaks
    ) +
    scale_y_discrete(
      name = "Episode length (days)",
      breaks = y_breaks, 
      labels = y_breaks_labs
    ) +
    scale_fill_distiller(
      name = "n patients\n(log10 scale)",
      trans = "log10",
      palette = "Blues",
      direction = 1
    ) +
    guides(
      fill = guide_colorbar(
        label =
      )
    ) +
    labs(
      subtitle = glue("Episode trigger event: {trigger}")
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90),
      axis.title.x = element_text(margin = margin(t=10)),
      axis.title.y = element_text(margin = margin(r=10)),
      legend.position = c(0.75,0.75)
    )
  
  ggsave(
    filename = here::here("output", "eda", glue("episode_lengths_{trigger}.png")),
    width = 18, height = 16, units = "cm")
  
}

episode_length_plot()
for (i in c("covid_primary_care_code_date", "covid_primary_care_positive_test_date", "covidadmitted_date", "postest_date")) {
  episode_length_plot(i)
}


