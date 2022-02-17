################################################################################
# This script:
# - combines estimates from all models into csv for release

################################################################################
library(tidyverse)
library(RColorBrewer)
library(glue)

################################################################################
fs::dir_create(here::here("output", "report", "data"))

################################################################################
# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# define comparisons
comparisons <- c("BNT162b2", "ChAdOx", "both")

# redaction functions
source(here::here("analysis", "lib", "redaction_functions.R"))

################################################################################
model_tidy_list <- unlist(lapply(
  comparisons,
  function(x)
    unlist(lapply(
      subgroup_labels,
      function(y)
        lapply(
          unname(outcomes),
          function(z)
            try(
              readr::read_rds(
                here::here("output", "models_cox", "data", glue("modelcox_tidy_{x}_{y}_{z}.rds")
                )
              ) %>%
                mutate(comparison = x, subgroup = y, outcome = z)
            )
        )
    ),
    recursive = FALSE
    )
),
recursive = FALSE
)

model_tidy_tibble <- bind_rows(
  model_tidy_list[sapply(model_tidy_list, function(x) is_tibble(x))]
) %>%
  select(subgroup, comparison, outcome, model, variable, label, reference_row,
         n_obs, n_event,
         estimate, conf.low, conf.high) %>%
  mutate(across(model, 
                factor, levels = 1:2, labels = "unadjusted", "adjusted")) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) %>%
  # redact estimates where n_obs <= 5
  mutate(redact = redactor(n = n_obs, threshold = 5)) %>%
  mutate(across(c(n_obs, estimate, conf.low, conf.high), 
                ~if_else(redact, NA_real_, .x))) %>%
  # redact n_event where n_event <=5
  mutate(across(n_event, redactor2)) 
  
cat("\nEstimates redacted for the following:\n")
model_tidy_tibble %>%
  filter(redact) %>%
  distinct(subgroup, comparison, outcome, model, variable, label) %>%
  print(n=Inf)

readr::write_csv(
  model_tidy_tibble,
  here::here("output", "models_cox", "data", "estimates_all.csv"))
