################################################################################
# process comparisons data
library(tidyverse)
library(lubridate)
library(glue)

################################################################################
## import study_parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# import variable names
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds")
)

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))

################################################################################
# individuals eligible based on box c, d & e criteria 
# arm and split info
data_arm <- bind_rows(
  readr::read_rds(
    here::here("output", "data", "data_eligible_e_vax.rds")) %>%
    rename(arm=brand),
  readr::read_rds(
    here::here("output", "data", "data_eligible_e_unvax.rds")) %>%
    mutate(arm = "unvax")
)  %>%
  select(patient_id, arm, split)

# vars from data_processed
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) %>%
  select(patient_id, subgroup,
         jcvi_group, elig_date, region, 
         dereg_date, death_date, coviddeath_date, noncoviddeath_date,
         any_of(unname(model_varlist$demographic)))

# vax data
data_wide_vax_dates <- readRDS(
  here::here("output", "data", "data_wide_vax_dates.rds")) %>%
  select(patient_id, covid_vax_1_date, covid_vax_3_date)

# read data for ever covariates
data_covariates <- arrow::read_feather(
  file = here::here("output", "input_covs.feather")) 

################################################################################
# covid infection episodes

# covid events within 90 days of each other grouped into one covid episode
# BTW it may be useful to reduce this threshold (e.g. to 30)  when verifying
# that the below code works as desired
episode_length <- 90

### REVIEW UP TO LINE 130
data_episodes0 <- data_covariates %>%
  select(
    patient_id, 
    # select recurring events
    starts_with(c("postest", "covidadmitted", "covid_primary_care"))
    ) %>%
  # clean dates variables
  mutate(across(contains("_date"), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days"))) %>%
  # join covid death data
  left_join(data_processed %>% 
              select(patient_id, coviddeath_1_date = coviddeath_date),
            by = "patient_id") %>%
  # reshape
  pivot_longer(
    cols = -patient_id, 
    values_drop_na = TRUE
    ) %>%
  mutate(across(name, ~str_remove(.x, "_\\d+_date"))) %>%
  arrange(patient_id, value) %>%
  # by patient_id, calculate lag of date
  group_by(patient_id) %>%
  mutate(value_lag = lag(value)) %>%
  ungroup() %>%
  # calculate difference in date between subsequent events
  mutate(diff = as.integer(value - value_lag))

data_episodes_2plus <- data_episodes0 %>%
  # only keep events that occurred > episode_length days after previous event
  filter(diff > episode_length) %>%
  select(patient_id, name, value) %>%
  arrange(patient_id, value) %>%
  group_by(patient_id) %>%
  # assign episode number, starting at 2
  mutate(episode = row_number() + 1) %>%
  ungroup()

data_episodes <- data_episodes0 %>%
  left_join(data_episodes_2plus) %>%
  arrange(patient_id, value) %>%
  group_by(patient_id) %>%
  # assign episode number 1 to all events not in data_episodes_2plus
  # then calculate the cumultive maximum by patient_id
  mutate(across(episode, ~cummax(if_else(is.na(.x), 1, .x)))) %>%
  ungroup() %>%
  select(-diff, -value_lag) %>%
  # start and end date of each episode
  group_by(patient_id, episode) %>%
  mutate(
    episode_start_date = min(value),
    episode_end_date = max(value),
  ) %>%
  ungroup() %>%
  # first date of each event type in each episode
  group_by(patient_id, episode, name) %>%
  mutate(type_date = min(value)) %>%
  ungroup() %>%
  select(-value) %>%
  distinct() %>%
  # rehape
  pivot_wider(
    names_from = name,
    values_from = type_date,
    names_glue = "{name}_date"
  )

readr::write_rds(
  data_episodes,
  here::here("output", "data", "data_episodes.rds"),
  compress = "gz"
)

################################################################################
data_all <- data_arm %>%
  # join to covariates data
  left_join(
    data_covariates %>%
      select(patient_id, 
             matches(c("start_\\d_date", "end_\\d_date")),
             starts_with("anytest"), asplenia,
             any_of(unname(unlist(model_varlist)))) %>%
      mutate(across(contains("_date"), 
                    ~ floor_date(
                      as.Date(.x, format="%Y-%m-%d"),
                      unit = "days"))),
    by = "patient_id") %>%
  # join to data_processed
  left_join(
    data_processed, by = "patient_id"
  ) %>%
  # join to vaccines
  left_join(
    data_wide_vax_dates, 
    by = "patient_id"
  ) %>%
  # derive remaining covariates
  mutate(
    
    pregnancy = pregnancy & (sex == "Female") & (age < 50),
    
    immunosuppressed = immunosuppressed | asplenia,
    
    multimorb =
      # as.integer(bmi %in% "Obese III (40+)") +
      as.integer(chd)  +
      as.integer(diabetes) +
      as.integer(cld) +
      as.integer(ckd) +
      as.integer(crd) +
      as.integer(immunosuppressed) +
      as.integer(cns),
    
    multimorb = cut(
      multimorb,
      breaks = c(0, 1, 2, Inf),
      labels=c("0", "1", "2+"),
      right=FALSE)
    
  ) %>%
  mutate(across(test_hist_n,
                ~ factor(case_when(
                  is.na(.x) ~ NA_character_,
                  .x < 1 ~ "0",
                  .x < 2 ~ "1",
                  .x < 3 ~ "2",
                  TRUE ~ "3+"
                )))) %>%
  mutate(subsequent_vax_date = if_else(
    arm %in% "unvax",
    covid_vax_1_date,
    covid_vax_3_date)) %>%
  select(-covid_vax_1_date, -covid_vax_3_date, -asplenia)

readr::write_rds(
  data_all,
  here::here("output", "data", "data_all.rds"),
  compress = "gz"
)

################################################################################
# store min and max fu dates for each subgroup

# create output directory
fs::dir_create(here::here("output", "lib"))

# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

data_min_max_fu <- data_all %>%
  group_by(subgroup) %>%
  summarise(
    min_fu_date = min(start_1_date),
    max_fu_date = max(end_6_date),
    # round total to nereast 7 for disclosure control
    # n = ceiling_any(n(), to=7),
    .groups = "keep"
  ) %>% 
  ungroup() %>%
  mutate(across(max_fu_date,
                ~ pmin(as.Date(study_parameters$end_date), .x)))


# data for release
readr::write_csv(
  data_min_max_fu,
  here::here("output", "lib", "data_min_max_fu.csv")
)

