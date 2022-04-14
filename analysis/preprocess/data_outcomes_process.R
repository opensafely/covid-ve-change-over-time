library(tidyverse)
library(lubridate)

################################################################################
# ever data for outcomes
data_ever <- arrow::read_feather(
  file = here::here("output", "input_ever.feather")) %>%
  mutate(across(contains("_date"), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days"))) 
 
## create folders for outputs
fs::dir_create(here::here("output", "data"))

################################################################################
# derive hospitalisation variable
data_covidadmitted_new <- data_ever %>%
  select(patient_id, covidadmitted_post_date, 
         matches("covidadmitted_postest_\\d+_date")) %>%
  filter(!is.na(covidadmitted_post_date)) %>%
  pivot_longer(
    cols = starts_with("covidadmitted_postest"),
    names_pattern = "covidadmitted_postest_(.)_date",
    names_transform = as.integer,
    values_drop_na = FALSE
  ) %>%
  # for those with no positive tests,
  mutate(across(value, ~if_else(is.na(.x), covidadmitted_post_date, .x))) %>%
  # calculate days between each positive test and hospitalisation
  mutate(diff = as.numeric(value - covidadmitted_post_date)) %>%
  mutate(
    postest = case_when(
      is.na(diff) ~ "none",
      diff < -4*7 ~ ">4 weeks before",
      diff  > 0 ~ "after",
      TRUE ~ "0-4 weeks before"
    )) %>%
  mutate(across(value,
                ~case_when(
                  # if no positive test, impute hospitalisation date
                  postest == "none" ~ covidadmitted_post_date,
                  # if positive test >4 weeks before hospitalisation, impute hospitalisation date
                  postest == ">4 weeks before" ~ covidadmitted_post_date,
                  # if positive test after hospitalisation, impute hospitalisation date
                  postest == "after" ~ covidadmitted_post_date,
                  # otherwise unchanged
                  TRUE ~ .x
                ))) %>%
  # recalculate diff
  mutate(diff = as.numeric(value - covidadmitted_post_date)) %>%
  # for each patient, keep only the earliest date in the value column
  group_by(patient_id) %>%
  mutate(min_diff = min(diff)) %>%
  ungroup() %>%
  filter(diff == min_diff) %>%
  select(patient_id, covidadmitted_new_date = value, diff, postest)

################################################################################
# summarise the derivation methods
capture.output(
  data_covidadmitted_new %>%
    group_by(postest) %>%
    count() %>%
    ungroup() %>%
    mutate(percent = round(100*n/sum(n), 2)) %>%
    kableExtra::kable(),
  file = here::here("output", "eda", "counts_postest_covidadmitted.txt"),
  append = FALSE
)


################################################################################
# check distribution of postest when used
data_covidadmitted_new %>%
  filter(postest == "0-4 weeks before") %>%
  ggplot(aes(x = diff)) +
  stat_count() +
  labs(title = "Distribution of positive test dates when used for hospitalisation date") +
  theme_bw()
ggsave(
  filename = here::here("output", "eda", "dist_postest_covidadmitted.png")
)  

################################################################################
# derive final dataset and save
data_outcomes <- data_ever %>%
  select(patient_id, postest_date = positive_test_post_date) %>%
  full_join(
    data_covidadmitted_new %>%
      select(patient_id, covidadmitted_date = covidadmitted_new_date),
    by = "patient_id")

readr::write_rds(
  data_outcomes,
  here::here("output", "data", "data_outcomes.rds"),
  compress = "gz"
)


