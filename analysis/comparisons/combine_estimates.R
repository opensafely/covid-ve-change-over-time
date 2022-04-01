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
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds")
)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# define comparisons
comparisons <- c("BNT162b2", "ChAdOx1", "both")

# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
model_tidy_list <- unlist(lapply(
  comparisons,
  function(x)
    unlist(lapply(
      subgroup_labels,
      function(y)
        unlist(lapply(
          unname(outcomes),
          function(z)
            lapply(
              1:K,
              function(kk)
                try(
                  readr::read_rds(
                    here::here("output", "models_cox", "data", glue("modelcox_tidy_{x}_{y}_{z}_{kk}.rds")
                    )
                  ) %>%
                    mutate(comparison = x, subgroup = y, outcome = z, period = kk)
                )
            )
        ),
        recursive = FALSE
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
  select(subgroup, comparison, outcome, model, period, variable, label, reference_row,
         n_obs, n_event,
         estimate, conf.low, conf.high) %>%
  mutate(across(model, 
                factor, levels = 1:2, labels = c("unadjusted", "adjusted")))

readr::write_csv(
  model_tidy_tibble,
  here::here("output", "models_cox", "data", "estimates_all.csv"))
