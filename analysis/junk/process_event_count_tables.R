################################################################################

library(tidyverse)
library(glue)

################################################################################

# read table_events
table_events <- bind_rows(
  readr::read_csv(
    here::here("output", "tte", "tables", glue("event_counts_BNT162b2.csv"))),
  readr::read_csv(
    here::here("output", "tte", "tables", glue("event_counts_ChAdOx.csv"))) %>%
    filter(!(arm %in% "unvax"))
) %>%
  group_split(subgroup, outcome)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroups <- c(subgroups, "all")

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))
outcomes <- unname(outcomes)

if (file.exists(here::here("output", "tte", "tables", glue("tidy_events.txt")))) {
  file.remove(here::here("output", "tte", "tables", glue("tidy_events.txt")))
}

################################################################################
process_tables <- function(
  data
) {
  
  # define arms
  arms <- sort(unique(data$arm))
  # define subgroup
  subgroup <- unique(data$subgroup)
  subgroup_label <- which(subgroups == subgroup)
  # define outcome
  outcome <- unique(data$outcome)
  
  caption_string <- glue("Subgroup = {subgroup}; Outcome = {outcome}")
  
  # define column order
  cols_order <- lapply(
    arms,
    function(x)
      sapply(c("n", "person_years", "events"),
             function(y)
               glue("{x} {y}"))
  )
  cols_order <- unname(unlist(cols_order))
  
  # generate table
  table <- data %>% 
    select(-subgroup, -outcome) %>%
    distinct() %>%
    pivot_wider(
      names_from = arm,
      values_from = c(n, person_years, events),
      names_glue = "{arm} {.value}"
    ) %>%
    select(k, all_of(cols_order))
  
  # save
  capture.output(
    kableExtra::kable(
      table,
      format = "pipe",
      caption = caption_string
    ),
    file = here::here("output", "tte", "tables", glue("tidy_events.txt")),
    append=TRUE
  )
  
  
}

################################################################################

lapply(
  seq_along(table_events),
  function(x)
    process_tables(table_events[[x]])
)


