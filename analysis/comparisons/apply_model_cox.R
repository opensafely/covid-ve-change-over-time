################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "BNT162b2"
  outcome <- "postest"
  
} else{
  comparison <- args[[1]]
  outcome <- args[[2]]
}

################################################################################
# read data

fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "models"))

source(here::here("analysis", "lib", "model_cox.R"))

model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)

################################################################################

if (comparison %in% c("BNT162b2", "ChAdOx")) {
  
  trt <- str_c(comparison, "_vax")
  
  data <- readr::read_rds(
    here::here("output", "data", glue("data_tte_{comparison}_{outcome}.rds"))) %>%
    left_join(
      bind_rows(
        lapply(
          c("vax", "unvax"),
          function(x)
            readr::read_rds(
              here::here("output", "data", glue("data_comparisons_{comparison}_{x}.rds"))) %>%
            select(patient_id, comparison, jcvi_group, elig_date, region, 
                   unname(unlist(model_varlist)))
        )
      ),
      by = c("patient_id", "comparison"))
    
        
} else if (comparison %in% "both") {
  
  trt <- "BNT162b2_vax"
  
  data <- readr::read_rds(
    here::here("output", "data", glue("data_tte_{comparison}_{outcome}.rds"))) %>%
    left_join(
      bind_rows(
        lapply(
          c("BNT162b2", "ChAdOx"),
          function(x)
            readr::read_rds(
              here::here("output", "data", glue("data_comparisons_{x}_vax.rds"))) %>%
            select(patient_id, comparison, jcvi_group, elig_date, region, 
                   unname(unlist(model_varlist)))
        )
      ),
      by = c("patient_id", "comparison"))
  
}

################################################################################

data_strata <- data %>% 
  mutate(strata_var = factor(str_c(jcvi_group, elig_date, region, sep = ", ")))



################################################################################

events_per_personyears <- function(strata) {
  data_strata %>%
    filter(strata_var %in% strata) %>%
    mutate(days = tstop-tstart) %>%
    group_by(comparison, arm) %>%
    summarise(
      n = n(),
      personyears = sum(days)/365,
      events = sum(status),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    mutate(across(c(n, events, personyears),
                  # round to nearest 10
                  ~ scales::comma(round(.x, -1), accuracy = 1))) %>%
    mutate(value = str_c(events, " / ", personyears, " (",n,")")) %>%
    select(comparison, arm, value) %>%
    pivot_wider(names_from = arm, values_from = value) %>%
    kableExtra::kable(
      "pipe",
      caption = glue("Strata: {strata}; {outcome} events / person-years (n)")
    )
}

################################################################################

for (s in levels(data_strata$strata_var)) {
  print(events_per_personyears(s))
}

capture.output(
  for (s in levels(data_strata$strata_var)) {
    print(events_per_personyears(s))
  },
  file = here::here("output", "tables", glue("{comparison}_{outcome}_incidence.txt")),
  append=FALSE
)

################################################################################

data_cox <- data_strata %>%
  dummy_cols(
    select_columns = c("comparison"),
    remove_selected_columns = TRUE
  ) %>%
  mutate(across(starts_with("comparison"),
                ~ if_else(arm %in% trt,
                          .x, 0L))) %>%
  mutate(strata_var = factor(str_c(jcvi_group, elig_date, region, sep = ", ")))

################################################################################

model_output <- list()
model_output[[1]] <- cox_model(
  number = 0, 
  formula = formula_cox_0,
  filename_prefix = glue("{comparison}_{outcome}"))
model_output[[2]] <- cox_model(
  number = 1, 
  formula = formula_cox_1,
  filename_prefix = glue("{comparison}_{outcome}"))
model_output[[3]] <- cox_model(
  number = 2, 
  formula = formula_cox_2,
  filename_prefix = glue("{comparison}_{outcome}"))


model_summary <- bind_rows(
  lapply(
    # only bind tibbles (to avoid errors in case some models did not converge)
    seq_along(model_output)[sapply(model_output, function(x) is_tibble(x[[1]]))],
    # select summary
    function(x) model_output[[x]]$summary
  )) %>%
  mutate(outcome = outcome)
readr::write_rds(
  model_summary,
  here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_summary.rds"))) 

### postprocessing using broom (may be unreliable)
# combine results
model_glance <- bind_rows(
  lapply(
    # only bind tibbles (to avoid errors in case some models did not converge)
    seq_along(model_output)[sapply(model_output, function(x) is_tibble(x[[1]]))],
    # select glance
    function(x) model_output[[x]]$glance
  )) %>%
  mutate(outcome = outcome)

readr::write_csv(
  model_glance,
  here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_glance.csv"))) 

