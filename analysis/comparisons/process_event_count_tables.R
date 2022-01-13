################################################################################

library(tidyverse)
library(glue)

################################################################################
## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "BNT162b2"
  
} else{
  comparison <- args[[1]]
}

################################################################################

# read table_events
table_events <- readr::read_rds(
  here::here("output", "tte", "tables", glue("event_counts_{comparison}.rds")))

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroups <- c(subgroups, "all")

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))

################################################################################

process_tables <- function(
  .data
) {
  
  # bind across outcomes
  table <- bind_rows(.data)
  
  # define subgroup
  subgroup <- unique(table$subgroup)
  subgroup_label <- which(subgroups == subgroup)
  
  # save data for generating tables for paper
  readr::write_csv(
    table,
    here::here("output", "tte", "tables", glue("events_{comparison}_{subgroup_label}.csv")))
  
  # create tidy output for data checking
  table_list <- table %>% 
    mutate(across(outcome, factor, levels = outcomes)) %>%
    select(-subgroup) %>%
    group_split(outcome) %>%
    as.list()
  
  # write tidy output to txt
  capture.output(
    lapply(
      table_list,
      function(x) {
        outcome <-  unique(x$outcome)
        x %>%
          select(-outcome) %>%
          kableExtra::kable(
            "pipe",
            caption = glue("subgroup = {subgroup}, outcome = {outcome}")
          )
      }
    ),
    file = here::here("output", "tte", "tables", glue("tidy_events_{comparison}_{subgroup_label}.txt")),
    append=FALSE
  )
  
}

################################################################################

lapply(
  seq_along(table_events),
  function(x)
    process_tables(table_events[[x]])
)


