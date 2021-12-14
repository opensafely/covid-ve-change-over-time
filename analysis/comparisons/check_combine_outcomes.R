library(tidyverse)
library(glue)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobs <- FALSE
  group <- "02"
  
} else{
  removeobs <- TRUE
  group <- args[[1]]
}

# read data
data_outcomes <- readr::read_rds(
  here::here("output", glue("jcvi_group_{group}"), "data", "data_outcomes.rds")) 

# distribution of days between outcome events
plot_check <- data_outcomes %>%
  transmute(
    postest_covidadmitted = as.integer(covidadmitted_date - postest_date),
    postest_coviddeath = as.integer(coviddeath_date - postest_date),
    covidadmitted_coviddeath = as.integer(coviddeath_date - covidadmitted_date),
  ) %>%
  pivot_longer(
    cols = everything(),
    values_drop_na = TRUE
  ) %>%
  mutate(across(name,
         ~str_replace(.x, "_", " and "))) %>%
  group_by(name) %>%
  mutate(mean = mean(value)) %>%
  ungroup() %>%
  ggplot(aes(x = value, colour = name)) +
  geom_freqpoly() +
  geom_vline(aes(xintercept = mean, colour = name), linetype = "dashed") +
  labs(x = "days between events",
       caption = "Means are represented by dashed vertical lines.") +
  scale_color_discrete(name = "events") +
  coord_cartesian(ylim = c(5, NA))

ggsave(plot_check,
       filename = here::here("output", glue("jcvi_group_{group}"), "images", "check_combine_outcomes.png"),
       width=14, height=12, units="cm")

# how many have each combination of covid outcomes?
# note that this is across all comparisons
data_check <- data_outcomes %>%
  select(patient_id, postest_date, covidadmitted_date, coviddeath_date) %>%
  mutate(across(-patient_id, ~!is.na(.x))) %>%
  group_by(postest_date, covidadmitted_date, coviddeath_date) %>% 
  count() %>%
  ungroup() %>%
  # round n to closest 10
  mutate(n = round(n, -1))

readr::write_csv(
  data_check,
  here::here("output", glue("jcvi_group_{group}"), "tables", "check_combine_outcomes.csv"))
  