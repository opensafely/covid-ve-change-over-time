################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
library(gt)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "both"
  subgroup_label <- 5
  outcome <- "noncoviddeath"
  
} else{
  comparison <- args[[1]]
  subgroup_label <- args[[2]]
  outcome <- args[[3]]
}

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroups <- c(subgroups, "all")
subgroup <- subgroups[subgroup_label]

################################################################################
# read data

fs::dir_create(here::here("output", "models_cox", "data"))
fs::dir_create(here::here("output", "models_cox", "tables"))

model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)
vars <- unname(unlist(model_varlist))

# redaction functions
source(here::here("analysis", "lib", "redaction_functions.R"))

################################################################################
arm1 <- if_else(comparison == "ChAdOx", "ChAdOx", "BNT162b2")
arm2 <- if_else(comparison == "both", "ChAdOx", "unvax")

data_cox <- readr::read_rds(
  here::here("output", "tte", "data", glue("data_tte_{comparison}_{subgroup_label}_{outcome}.rds"))) %>%
  left_join(
    bind_rows(
      readr::read_rds(
        here::here("output", "comparisons", "data", glue("data_comparisons_{arm1}.rds"))),
      readr::read_rds(
        here::here("output", "comparisons", "data", glue("data_comparisons_{arm2}.rds")))
    ) %>%
      select(patient_id, comparison, jcvi_group, elig_date, region, 
             unname(unlist(model_varlist))),
    by = c("patient_id", "comparison")) %>% 
  mutate(strata_var = factor(str_c(jcvi_group, elig_date, region, sep = ", "))) %>%
  droplevels()


readr::write_rds(
  data_cox,
  here::here("output", "models_cox", "data", glue("data_cox_{comparison}_{subgroup_label}_{outcome}.rds")),
  compress = "gz"
)

################################################################################

total_events <- data_cox %>% filter(status) %>% nrow()
model_instructions <- list(
  model = total_events > 0
)

readr::write_rds(
  model_instructions,
  here::here("output", "lib", glue("model_instructions_{comparison}_{subgroup_label}_{outcome}.rds"))
)

if (model_instructions$model) {
  
  cat("...split data by comparison and status...\n")
  tbl_list <- data_cox %>%
    select(comparison, status, all_of(vars)) %>%
    group_split(comparison, status) 
  
  # names for each element in list
  group_split_labels <- lapply(
    tbl_list,
    function(x) str_c(unique(x$comparison), unique(x$status))
  ) %>% unlist()
  
  cat("...summarise number of events...\n")
  # summarise the number of events by level of covariates (within comparisons)
  tbltab_list <- tbl_list %>%
    map(~.[,-c(1,2)]) %>%
    map(
      function(data){
        map(data,
            function(x) {
              tab <- table(x)
              tibble(.level = names(tab), 
                     n = as.vector(tab))
            }) %>%
          bind_rows(.id="variable") 
      }
    )
  
  # apply names
  names(tbltab_list) <- group_split_labels
  
  cat("...prepare table...\n")
  tbltab <- bind_rows(
    tbltab_list,
    .id = "group"
  ) %>%
    mutate(
      comparison = str_extract(group, "\\d"),
      status = as.logical(str_remove(group, "\\d")))  %>%
    select(-group) %>%
    pivot_wider(
      names_from = status,
      values_from = n
    ) %>%
    pivot_wider(
      names_from = comparison,
      values_from = c("FALSE", "TRUE"),
      names_glue = "comparison{comparison}_{.value}"
    )
  
  cat("...format and save table...\n")
  tbltab %>%
    gt(
      groupname_col="variable",
      rowname_col = ".level"
    ) %>%
    tab_spanner_delim("_") %>%
    tab_stubhead(label = "variable") %>%
    opt_css(css = ".gt_stub { padding-left: 50px !important; }") %>%
    fmt_number(
      columns = starts_with(c("comparison")),
      sep_mark = ",",
      decimals = 0
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = "lightcyan")
      ),
      locations = cells_body(
        columns = ends_with("TRUE")
      )
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = "lightcyan")
      ),
      locations = cells_column_labels(
        columns = ends_with("TRUE")
      )
    ) %>%
    gtsave(
      filename = glue("eventcheck_{comparison}_{subgroup_label}_{outcome}.html"),
      path = here::here("output", "models_cox", "tables")
    )
  
} else {
  
  readr::write_file(
    x="",
    here::here("output", "models_cox", "tables", glue("eventcheck_{comparison}_{subgroup_label}_{outcome}.html")),
    append = FALSE
    )
  
}


 