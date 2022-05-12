################################################################################
# check distribution of number of positive tests, probably covid recordings and hospitalisations
# this will be used to determine how many recurring variables to define in study_definition_covs

library(tidyverse)

# create output directory
fs::dir_create(here::here("output", "eda"))

# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

data_processed %>%
  select(patient_id, ends_with("_n")) %>%
  pivot_longer(cols = -patient_id) %>%
  group_by(name, value) %>%
  count() %>%
  ungroup() %>%
  mutate(across(n, ~ceiling_any(.x, to = 7))) %>%
  # filter(value < 50) %>% # to avoid very large numbers messing up the scales.
  ggplot(aes(x = value, y = n)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~name, nrow=5, scales = "free")

ggsave(
  filename = here::here("output", "eda", "prior_covid_outcomes_n.png"),
  width = 14, height = 20, units = "cm"
)
