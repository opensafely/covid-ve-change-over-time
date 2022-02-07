################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
library(fastDummies)
library(gt)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "both"
  subgroup_label <- 3
  outcome <- "postest"
  
} else{
  comparison <- args[[1]]
  subgroup_label <- args[[2]]
  outcome <- args[[3]]
}

################################################################################
# create output directories
fs::dir_create(here::here("output", "preflight", "data"))
fs::dir_create(here::here("output", "preflight", "tables"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroups <- c(subgroups, "all")
subgroup <- subgroups[subgroup_label]

# model covariates
model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds")
)
vars <- unname(unlist(model_varlist))

################################################################################
# read functions

# redaction functions
source(here::here("analysis", "lib", "redaction_functions.R"))

################################################################################
arm1 <- if_else(comparison == "ChAdOx", "ChAdOx", "BNT162b2")
arm2 <- if_else(comparison == "both", "ChAdOx", "unvax")

data_0 <- readr::read_rds(
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

################################################################################
# tabulate events per level

# check there are > 0 events
total_events <- data_0 %>% filter(status) %>% nrow()

if (total_events > 0) {
  
  cat("...split data by comparison and status...\n")
  tbl_list <- data_0 %>%
    select(comparison, status, all_of(vars)) %>%
    select(-age) %>%
    group_split(comparison, status)
  
  # names for each element in list
  group_split_labels <- lapply(
    tbl_list,
    function(x) str_c(unique(x$comparison), unique(x$status))
  ) %>% 
    unlist()
  
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
    ) %>%
    mutate(across(starts_with("comparison"), 
                  ~ if_else(is.na(.x), 0L, .x))) %>% 
    group_by(variable) %>%
    mutate(across(starts_with("comparison"), redactor2)) %>% 
    mutate(across(starts_with("comparison"), 
                  ~ if_else(is.na(.x), "-", scales::comma(.x, accuracy = 1)))) %>% 
    ungroup()
  
  cat("...format and save table...\n")
  tbltab %>%
    gt(
      groupname_col="variable",
      rowname_col = ".level"
    ) %>%
    tab_spanner_delim("_") %>%
    tab_stubhead(label = "variable") %>%
    opt_css(css = ".gt_stub { padding-left: 50px !important; }") %>%
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
      filename = glue("eventcheck_{comparison}_{subgroup_label}_{outcome}_REDACTED.html"),
      path = here::here("output", "preflight", "tables")
    )
  
  ################################################################################
  # remove comparisons with <= 10 events
  events_threshold <- 10
  
  # check events per comparison
  events_per_comparison <- data_0 %>%
    group_by(comparison) %>%
    summarise(events = sum(status), .groups="keep") %>%
    ungroup() %>%
    mutate(keep = events > events_threshold)
  
  keep_comparisons <- as.integer(events_per_comparison$comparison[events_per_comparison$keep])
  drop_comparisons <- as.integer(events_per_comparison$comparison[!events_per_comparison$keep])
  
  data_1 <- data_0 %>%
    filter(comparison %in% keep_comparisons) %>%
    droplevels()
  
  ################################################################################
  # check levels per covariate
  n_levels <- sapply(
    data_1 %>% select(all_of(vars)) %>% mutate(across(everything(), as.factor)),
    function(x) length(levels(x))
  )
  
  # calculate events per level
  events_per_level <- data_1 %>%
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
      # remove if one level or binary and one level with too few events
      n_levels ==1 | (n_levels == 2  & n_keep == 1) 
    ) %>% 
    distinct(variable)
  
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
          mutate(
            # number of events per new level
            new_n = sum(n),
            # index for new level
            min_index = min(index)
          ) %>%
          ungroup() %>%
          # update keep
          mutate(keep = new_n > threshold)
        
        # first min_index with <= 5 events
        first_false <- min(data_old$min_index[!data_old$keep])
        # unique values of min_index
        unique_min_index <- unique(data_old$min_index)
        # merge up unless first_false is the top level, in which case merge down
        if (first_false < max(unique_min_index)) {
          # merge up
          merge_levs <- c(first_false, unique_min_index[c(which(unique_min_index == first_false)+1)])
        } else {
          # merge down
          merge_levs <- c(unique_min_index[c(which(unique_min_index == first_false)-1)], first_false)
        }
        # merge the labels
        merged_lev <- str_c(data_old$new_level[merge_levs], collapse = " & ")
        
        # merge levels and add labels
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
  
  # apply the merge function
  new_level_key <- oridinal_var_list %>%
    map(~merge_levs_fun(., threshold = 10)) 
  
  # use the new merged levels to re-code the variables in the original data
  data_2 <- data_1
  for (i in seq_along(new_level_key)) {
    
    if (!is_empty(new_level_key[[i]])) {
      
      var <- unique(new_level_key[[i]]$variable)
      levs <- unique(new_level_key[[i]]$new_level)
      
      join_by <- "level"
      names(join_by) <- var
      
      data_2 <- data_2 %>%
        mutate(across(all_of(var), as.character)) %>%
        left_join(
          new_level_key[[i]] %>% select(-variable),
          by = join_by
        ) %>%
        mutate(!! sym(var) := factor(new_level, levels = levs)) %>%
        select(-new_level)
      
    } 
    
  }
  
  # if merge resulted in 1 level, drop the variable
  drop_merged_var <- sapply(
    new_level_key,
    function(x) {
      if (is_empty(x)) {
        return(NA_character_)
      } else {
        drop <- n_distinct(x$new_level) == 1
        if (drop) return(unique(x$variable)) else return(NA_character_)
      }
    }
  )
  
  drop_merged_var <- drop_merged_var[!is.na(drop_merged_var)]
  
  data_3 <- data_2 %>%
    select(-all_of(c(drop_vars$variable, drop_merged_var))) %>%
    droplevels()
  
  ################################################################################
  # create comparison dummy variables
  data_4 <- data_3 %>%
    dummy_cols(
      select_columns = c("comparison"),
      remove_selected_columns = TRUE
    ) %>%
    mutate(across(starts_with("comparison"),
                  ~ if_else(arm %in% arm1,
                            .x, 0L))) %>%
    rename_at(vars(starts_with("comparison")), ~str_c(.x, "_", arm1))
  
  ################################################################################
  # create age variables
  if (subgroup_label == 1) {
    
    # age and age^2 for subgroup 16-64 and vulnerable
    data_5 <- data_4 %>%
      mutate(
        age_1 = age,
        age_1_squared = age^2
        ) %>%
      select(-age)
    
  } else {
    
    data_4_list <- data_4 %>%
      group_split(jcvi_group, elig_date) %>%
      as.list()
    
    for (i in seq_along(data_4_list)) {
      
      data_4_list[[i]] <- data_4_list[[i]] %>%
        mutate(!! sym(glue("age_{i}")) := age) %>%
        select(-age)
      
      # add an age^2 term for jcvi group 2 (80+)
      if (unique(data_4_list[[i]]$jcvi_group) == "02") {
        
        data_4_list[[i]] <- data_4_list[[i]] %>%
          mutate(!! sym(glue("age_{i}_squared")) := sym(glue("age_{i}"))^2)
        
      }
      
    }
    
    data_5 <- bind_rows(
      data_4_list
    ) %>%
      mutate(across(starts_with("age"),
                    ~ if_else(is.na(.x),
                              0,
                              as.double(.x))))
    
  }
  
  ################################################################################
  # define formulas
  
  comparisons <- data_5 %>% select(starts_with("comparison")) %>% names()
  
  formula_unadj <- as.formula(str_c(
    "Surv(tstart, tstop, status, type = \"counting\") ~ ",
    str_c(comparisons, collapse = " + "),
    " + strata(strata_var)"))
  
  demog_vars <- c(names(data_5)[str_detect(names(data_5), "^age_")],
                  unname(model_varlist$demographic[which(model_varlist$demographic %in% names(data_5))]))
  formula_demog <- as.formula(str_c(c(". ~ . ", demog_vars), collapse = " + "))
  
  clinical_vars <- unname(model_varlist$clinical[which(model_varlist$clinical %in% names(data_5))])
  formula_clinical <- as.formula(str_c(c(". ~ . ", clinical_vars), collapse = " + "))
  
  ################################################################################
  
  formulas_list <-  list(
    "unadjusted" = formula_unadj, 
    "demographic" = formula_demog, 
    "clinical" = formula_clinical)
  
  model_input <- list(
    data = data_5,
    formulas = formulas_list
  )
  
  readr::write_rds(
    model_input,
    here::here("output", "preflight", "data", glue("model_input_{comparison}_{subgroup_label}_{outcome}.rds"))
  )
  
  ################################################################################
  
  preflight_report <- function(
    dropped_comparisons,
    dropped_variables,
    merged_variables,
    subgroup_string = subgroup
  ) {
    ####
    cat(glue("Comparison = {comparison}; Subgroup = {subgroup_string}; Outcome = {outcome}"), "\n")
    cat("---\n")
    if (is_empty(drop_comparisons)) {
      dropped_comparisons <- "none"
    } else {
      dropped_comparisons <- str_c(dropped_comparisons, collapse = ", ")
    }
    cat(glue("Dropped comparisons: {dropped_comparisons}"), "\n")
    ####
    cat("---\n")
    if (is_empty(dropped_variables)) {
      dropped_comparisons <- "none"
    } else {
      dropped_variables <- str_c(str_c("- ", dropped_variables), collapse = "\n")
    }
    cat(glue("Dropped variables:\n{dropped_variables}"), "\n")
    ####
    cat("---\n")
    if (is_empty(merged_variables)) {
      cat("No levels merged.", "\n")
    } else {
      
      merged_variables %>%
        kableExtra::kable("pipe",
                          caption = "Merged levels:") %>%
        print()
    }
    ####
    cat("\n")
    cat("---\n")
    cat(glue("Formulas:"), "\n")
    print(formulas_list)
    
  }
  
  capture.output(
    preflight_report(
      dropped_comparisons = drop_comparisons,
      dropped_variables = c(drop_vars$variable,drop_merged_var),
      merged_variables = bind_rows(new_level_key)
    ),
    file = here::here("output", "preflight", "tables", glue("preflight_report_{comparison}_{subgroup_label}_{outcome}.txt")),
    append = FALSE
  )
  
} else {
  # empty outputs to avid errors
  
  readr::write_file(
    x="",
    here::here("output", "preflight", "tables", glue("eventcheck_{comparison}_{subgroup_label}_{outcome}_EMPTY.html")),
    append = FALSE
  )
  
  readr::write_rds(
    tibble(),
    here::here("output", "preflight", "data", glue("model_input_{comparison}_{subgroup_label}_{outcome}.rds"))
  )
  
  capture.output(
    print("No events"),
    file = here::here("output", "preflight", "tables", glue("preflight_report_{comparison}_{subgroup_label}_{outcome}.txt")),
    append = FALSE
  )
  
}
