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
  here::here("output", "lib", "study_parameters.rds"))

outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# define comparisons
comparisons <- c("BNT162b2", "ChAdOx", "both")

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
                here::here("output", "models_cox", "data", glue("modelcox_summary_{x}_{y}_{z}.rds")
                )
              ) %>%
                mutate(comparison = x, subgroup = y)
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
  select(subgroup, comparison, outcome, model, term, everything())

readr::write_csv(
  model_tidy_tibble,
  here::here("output", "report", "data", "hr_all.csv"))

model_tidy_tibble_arms <- model_tidy_tibble %>%
    filter(
      str_detect(term, "^comparison")
      ) 

readr::write_csv(
  model_tidy_tibble_arms,
  here::here("output", "report", "data", "hr_arms.csv"))
