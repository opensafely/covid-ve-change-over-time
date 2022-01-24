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
  comparison <- "BNT162b2"
  subgroup_label <- 1
  outcome <- "postest"
  
} else{
  comparison <- args[[1]]
  subgroup_label <- args[[2]]
  outcome <- args[[3]]
}

# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

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

data_all <- readr::read_rds(
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


# readr::write_rds(
#   data_cox,
#   here::here("output", "models_cox", "data", glue("data_cox_{comparison}_{subgroup_label}_{outcome}.rds")),
#   compress = "gz"
# )

################################################################################
events_threshold <- 5

# get rid of comparisons with <= 5 events
events_per_comparison <- data_all %>%
  group_by(comparison) %>%
  summarise(events = sum(status)) %>%
  ungroup() %>%
  mutate(keep = events > events_threshold)

keep_comparisons <- as.integer(events_per_comparison$comparison[events_per_comparison$keep])
drop_comparisons <- as.integer(events_per_comparison$comparison[!events_per_comparison$keep])

data_comp <- data_all %>%
  filter(comparison %in% keep_comparisons) %>%
  droplevels()

# check levels 
n_levels <- sapply(
  data_comp %>% select(all_of(vars)) %>% mutate(across(everything(), as.factor)),
  function(x) length(levels(x))
)

# calculate events per level
events_per_level <- data_comp %>%
  filter(status) %>%
  select(all_of(vars)) %>%
  map(function(x) {
        tab <- table(x)
        tibble(level = names(tab), 
               n = as.vector(tab))
      }) %>%
  bind_rows(.id="variable") %>%
  mutate(keep = n > events_threshold) %>%
  group_by(variable) %>%
  mutate(
    n_levels = n(),
    n_keep = sum(keep)
    ) %>%
  ungroup()

drop_vars <- events_per_level %>%
  filter(
    # remove if one level or binary and one level with < threshold events
    n_levels ==1 | (n_levels == 2  & n_keep == 1) 
  ) %>% 
  select(variable)

# for ordinal variables, try combining levels 
oridinal_var_list <- events_per_level %>%
  filter(n_levels > 2) %>%
  group_split(variable)


# merge levels for ordinal variables with > 2 levels, in which some levels have low numbers
merge_levs_fun <- function(data, threshold=10) {
  
  # flag n that are less than threshold
  data <- data %>% mutate(keep = n > threshold)
  
  if (all(data$keep)) {
    
    # return empty tibble if all less that threshold
    return(tibble())
    
  } else {
    
    merge_levs_i <- function(data_in) {
      
      if (all(data_in$keep)) stop(glue("All levels have greater than {threshold} events."))
      
      # only applied in first loop
      if (!"new_level" %in% names(data_in)) {
        data_in <- data_in %>% 
          mutate(
            new_level = level,
            index = row_number()
          )
      }
      
      data_old <- data_in %>% 
        group_by(new_level) %>%
        mutate(new_n = sum(n),
               min_index = min(index)) %>%
        ungroup() %>%
        mutate(keep = new_n>threshold)
      
      # first index with <= 5 events
      first_false <- min(data_old$min_index[!data_old$keep])
      unique_min_index <- unique(data_old$min_index)
      # merge up unless last level, in which case merge down
      if (first_false < max(unique_min_index)) {
        merge_levs <- c(first_false,  unique_min_index[c(which(unique_min_index == first_false)+1)])
      } else {
        merge_levs <- c( unique_min_index[c(which(unique_min_index == first_false)-1)], first_false)
      }
      # merge the labels
      merged_lev <- str_c(data_old$new_level[merge_levs], collapse = " & ")
      
      # merge levels and apply labels
      data_new <- data_old %>%
        mutate(across(new_level, 
                      ~if_else(min_index %in% merge_levs, merged_lev, .x))) %>%
        group_by(new_level) %>%
        mutate(new_n = sum(n)) %>%
        ungroup() %>%
        mutate(keep = new_n>threshold) %>%
        select(-new_n)
      
      return(data_new)
      
    }
    
    for (i in 1:(length(data$level)-1)) { # max number of possible merges 
      
      if (!all(data$keep)) {
        data <- merge_levs_i(data)
      } 
      
    }
    
    return(data %>% select(variable, level, new_level))
    
  }
  
}



new_level_key <- oridinal_var_list %>%
  map(~merge_levs_fun(., threshold = 10)) 

# use the new merged levels to recode the variables in the original data
for (i in seq_along(new_level_key)) {
  
  if (!is_empty(new_level_key[[i]])) {
    
    var <- unique(new_level_key[[i]]$variable)
    levs <- unique(new_level_key[[i]]$new_level)
    
    join_by <- "level"
    names(join_by) <- var
    
    data_comp <- data_comp %>%
      mutate(across(all_of(var), as.character)) %>%
      left_join(
        new_level_key[[i]] %>% select(-variable),
        by = join_by
      ) %>%
      mutate(!! sym(var) := factor(new_level, levels = levs)) %>%
      select(-new_level)
    
  } 
  
}


data_comp %>%
  select(-all_of(drop_vars$variable))




################################################################################
# check there are >0 events
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


 