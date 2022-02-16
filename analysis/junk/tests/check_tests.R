# run checks on input_tests.feather

## setup
library(tidyverse)

## source functions
source(here::here("analysis", "lib", "data_properties.R"))

## create folders for outputs
fs::dir_create(here::here("output", "tests", "images"))
fs::dir_create(here::here("output", "tests", "tables"))


data_tests <- arrow::read_feather(
  file = here::here("output", "input_tests.feather")) %>%
  mutate(across(c(contains("_date")), 
                ~ floor_date(
                  as.Date(., format="%Y-%m-%d"),
                  unit = "days")))

# tabulate all vars
data_properties(
  data = data_tests,
  path = file.path("output", "tests", "tables")
)

# coviariates
data_tests %>%
  select(patient_id, 
         pre=covid_test_pre_elig_n,
         post=covid_test_post_elig_n) %>%
  mutate(both = pre + post) %>%
  pivot_longer(cols = -patient_id) %>%
  ggplot(aes(x = value)) +
  geom_bar() +
  facet_wrap(~ name, scales = "free") +
  theme(legend.position = "bottom")
ggsave(
  filename = here::here("output", "tests", "images", "covariate_distribution.png"),
  width=20, height=15, units="cm")

