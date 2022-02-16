################################################################################

library(tidyverse)
library(glue)

################################################################################

# read table_events
table_events <- bind_rows(
  readr::read_rds(
    here::here("output", "tte", "tables", glue("event_counts_BNT162b2.rds"))),
  readr::read_rds(
    here::here("output", "tte", "tables", glue("event_counts_ChAdOx.rds"))) %>%
    filter(!(arm %in% "unvax"))
) %>%
  group_split(subgroup)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroups <- c(subgroups, "all")

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))
outcomes <- unname(outcomes)

# redaction functions
source(here::here("analysis", "lib", "redaction_functions.R"))

################################################################################
process_tables <- function(
  data
) {
  
  # define subgroup
  subgroup <- unique(data$subgroup)
  subgroup_label <- which(subgroups == subgroup)
  
  caption_string <- glue("Subgroup = {subgroup}")
  
  # define column order
  cols_order <- lapply(
    c("BNT162b2", "ChAdOx", "unvax"),
    function(x)
      sapply(c("n", outcomes),
             function(y)
               glue("{x}_{y}"))
  )
  cols_order <- unname(unlist(cols_order))
  
  # generate table
  table <- data %>% 
    distinct() %>%
    pivot_longer(cols = c(n, events)) %>%
    mutate(col_names = glue("{arm}_{outcome}_{name}")) %>%
    select(-arm, -outcome, -name) %>%
    distinct() %>%
    pivot_wider(names_from = col_names, values_from = value) %>%
    rename_at(vars(matches(".+_noncoviddeath_n")), ~ str_remove(.x, "_noncoviddeath")) %>%
    select(-matches(".+_.+_n")) %>%
    rename_at(vars(ends_with("_events")), ~ str_remove(.x, "_events")) %>%
    select(-subgroup)  
  
  # for subgroups in which only one vaccine used
  cols_order <- cols_order[cols_order %in% names(table)]
  
  # process table col names
  col_names_tidy <- str_replace(cols_order, "_", " ")
  col_names_tidy <- str_replace(col_names_tidy, "anytest", "any test")
  col_names_tidy <- str_replace(col_names_tidy, "postest", "+ve test")
  col_names_tidy <- str_replace(col_names_tidy, "covidadmitted", "C-19 hosp. admission")
  col_names_tidy <- str_replace(col_names_tidy, "noncoviddeath", "Non-C-19 death")
  col_names_tidy <- str_replace(col_names_tidy, "coviddeath", "C-19 death")
  
  # redact low values
  table_out <- table %>%
    select(k, all_of(cols_order)) %>%
    mutate(across(-k, 
                  ~ if_else(.x <= 5,
                            "-",
                            scales::comma(.x, accuracy = 1)))) 
  
  # save
  capture.output(
    kableExtra::kable(
      table_out,
      format = "pipe",
      col.names = c("k", col_names_tidy),
      caption = caption_string
    ),
    file = here::here("output", "tte", "tables", glue("tidy_events_{subgroup_label}.txt")),
    append=FALSE
  )
  
  
}

################################################################################

lapply(
  seq_along(table_events),
  function(x)
    process_tables(table_events[[x]])
)


