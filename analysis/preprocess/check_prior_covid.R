################################################################################
# check distribution of number of positive tests, probably covid recordings and hospitalisations
# this will be used to determine how many recurring variables to define in study_definition_covs

library(tidyverse)

# create output directory
fs::dir_create(here::here("output", "eda"))
fs::dir_create(here::here("output", "lib"))

# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

data_plot <- data_processed %>%
  select(patient_id, ends_with("_n")) %>%
  pivot_longer(cols = -patient_id) %>%
  group_by(name, value) %>%
  count() %>%
  ungroup() %>%
  arrange(name, desc(value)) %>%
  group_by(name) %>%
  mutate(
    cml_n = cumsum(n),
    thresh_value = if_else(
      # choose threshold to keep all tests for at least 99.9% of participants
      cml_n > 0.001*nrow(data_processed),
      value,
      NA_integer_
      )
    ) %>%
  mutate(across(thresh_value, max, na.rm = TRUE))  %>%
  ungroup() 


data_plot %>%
  mutate(across(n, ~log(ceiling_any(.x, to = 7)))) %>%
  filter(value <=20) %>% # to avoid very large numbers messing up the scales.
  ggplot(aes(x = value, y = n)) +
  geom_bar(stat = "identity", width = 1) +
  geom_vline(aes(xintercept = thresh_value), colour = "magenta", linetype = "dashed") +
  facet_wrap(~name, nrow=5) +
  labs(x = "n tests (truncated at 20)", y = "log(n) patients") +
  theme_bw()

ggsave(
  filename = here::here("output", "eda", "prior_covid_outcomes_n.png"),
  width = 14, height = 20, units = "cm"
)

################################################################################
# save thresholds for use in study definition
data_thresholds <- data_plot %>%
  distinct(name, thresh_value)

prior_covid_n <- as.list(data_thresholds$thresh_value)
names(prior_covid_n) <- data_thresholds$name

jsonlite::write_json(
  prior_covid_n, 
  path = here::here("output", "lib", "prior_covid_n.json"), 
  auto_unbox = TRUE, pretty=TRUE)
