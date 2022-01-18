library(tidyverse)
library(glue)
library(gt)


# redaction functions
source(here::here("analysis", "lib", "redaction_functions.R"))

data <- data_cox

formula <- update(status ~ 1, formula_demog) %>% update(formula_clinical)

stratum <- "comparison"

septab <- function(data, formula, stratum, brand, outcome, name){
  
  tbltab <- data %>%
    select(stratum, all.vars(formula)) %>%
    mutate(
      across(
        where(is.integer),
        ~as.character(.)
      )
    ) %>%
    split(.[[1]]) %>%
    # map(~.[,-1] %>% select(stratum, everything())) %>%
    map(
      function(data){
        map(data, redacted_summary_cat, redaction_threshold=0) %>%
          bind_rows(.id="variable") %>%
          select(-redacted, -pct_nonmiss)
      }
    )
      
  tbltab %>%
    bind_rows(.id = "comparison") %>%
    pivot_wider(
      id_cols=c(variable, .level),
      names_from = comparison,
      names_glue = "comparison{comparison}_{.value}",
      values_from = c(n, pct)
    ) %>%
    select(variable, .level, starts_with("comparison")) %>%
    filter(!str_detect(variable, "^comparison")) %>%
    gt(
      groupname_col="variable",
    ) %>%
    tab_spanner_delim("_") %>%
    fmt_number(
      columns = ends_with(c("pct")),
      decimals = 1,
      scale_by=100,
      pattern = "({x})"
    ) #%>%
    gtsave(
      filename = glue("sepcheck_{stratum}_{recentpostest_period}_{brand}_{outcome}_{name}.html"),
      path=here("output", cohort, "descriptive", "model-checks", strata_var, recentpostest_period)
    )
}