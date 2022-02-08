library(tidyverse)
library(glue)
library(lubridate)
library(gt)

################################################################################

fs::dir_create(here::here("output", "report", "tables"))

################################################################################
## import study_parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))


# individuals eligible based on box c & e criteria
data_eligible_e_vax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_vax.rds"))  %>%
  mutate(arm = brand) %>%
  select(patient_id, start_of_period, arm)

# individuals eligible based on box d & e criteria
data_eligible_e_unvax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_unvax.rds"))  %>%
  mutate(arm = "unvax") %>%
  select(patient_id, start_of_period, arm)

# read list of covariates for model
model_varlist <- readr::read_rds(
  here::here("output", "lib", "model_varlist.rds"))

# read strata_vars
strata_vars <- readr::read_rds(
  here::here("output", "lib", "strata_vars.rds"))
strata_vars <- strata_vars[strata_vars!="elig_date"]

# processed data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) 

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))
outcomes <- unname(outcomes)
outcomes <- outcomes[outcomes!="anytest"]

# read script for processing covariates for comparison 1
source(here::here("analysis", "lib", "process_covariates.R"))

# redaction functions
source(here::here("analysis", "lib", "redaction_functions.R"))

################################################################################
# function to be applied in dplyr::filter
no_evidence_of <- function(cov_date, index_date) {
  is.na(cov_date) | index_date < cov_date
}

censor_vars <- c("death_date", "dereg_date")

data_comparison_1 <- bind_rows(
  data_eligible_e_vax,
  data_eligible_e_unvax
  ) %>%
  left_join(data_processed,
            by = "patient_id") %>%
  mutate(
    # start date of comparison 1 
    start_fu_date = start_of_period + days(14),
    end_fu_date = start_fu_date + days(28)
  ) %>%
  select(-start_of_period) %>%
  # remove if death or dereg before start_of_period
  filter_at(
    all_of(censor_vars),
    all_vars(no_evidence_of(., start_fu_date))) %>%
  ungroup() %>%
  select(patient_id, jcvi_group, elig_date, region, arm,
         start_fu_date, end_fu_date) %>%
  mutate(comparison = factor(1, levels = 1:study_parameters$max_comparisons))

################################################################################

data_tables <- data_comparison_1 %>%
  process_covariates() %>%
  select(patient_id, arm, region, jcvi_group, subgroup,
         all_of(unname(unlist(model_varlist)))) %>% 
  group_split(subgroup)

################################################################################
# make table1 for all and each subgroup
for (i in c(0, seq_along(data_tables))) {
  cat(glue("---- loop {i} ----"), "\n")
  cat("---- define obejcts ----\n")
  variables <- c(unname(strata_vars), unname(unlist(model_varlist)))
  variables <- variables[variables != "age"]
  vars_ordered_levs <- c("region", "jcvi_group", "sex", "imd", "ethnicity", "bmi", "multimorb", "test_hist_1_n")
  
  # tibble for assigning tidy variable names
  var_labels <- tibble(
    variable = c(strata_vars, model_varlist$clinical, model_varlist$demographic),
    variable_label = names(c(strata_vars, model_varlist$clinical, model_varlist$demographic))
  )
  
  if (i == 0) {
    
    data <- bind_rows(data_tables) %>%
      droplevels()
    
    subgroup <- "All subgroups"
    subgroup_label <- 5
    
    variables <- c("subgroup", variables)
    vars_ordered_levs <- c("subgroup", vars_ordered_levs)
    
    var_labels <- var_labels %>%
      add_row(variable = "subgroup", variable_label = "Subgroup",
              .before=TRUE)
    
    min_elig_date <- "2020-12-08"
    
  } else {
    
    data <- data_tables[[i]] %>%
      droplevels()
    
    subgroup <- unique(data$subgroup)
    subgroup_label <- which(subgroups == subgroup)
    
    min_elig_date <- data_processed %>%
      filter(subgroup %in% subgroup) %>%
      summarise(min_elig_date = min(elig_date))
    min_elig_date <- min_elig_date$min_elig_date
    
  }
  
  # function for creating tibble of categories for each variable
  var_tibble <- function(var) {
    var_region <- tibble(
      variable = var,
      category = levels(data[[var]])
    )
  }
  # tibble for specifying order of variables and categories
  var_order <- tibble(
    variable = variables
  ) %>%
    left_join(
      bind_rows(
        lapply(vars_ordered_levs,
               var_tibble)),
      by = "variable"
    ) %>%
    mutate(across(category, ~ if_else(is.na(.x), "yes", .x)))
  
  cat("---- summarise variables ----\n")
  # summarise each variable and (within variables) redact values <=5
  summary_var <- function(.data, var) {
    out <- .data %>%
      group_by(arm, !! sym(var)) %>%
      count() %>%
      ungroup(!! sym(var)) %>%
      mutate(arm_total = sum(n)) %>%
      ungroup() %>%
      mutate(percent = round(100*n/arm_total,0)) %>%
      group_by(arm, !! sym(var)) %>%
      mutate(across(n, redactor2)) %>%
      ungroup() %>%
      mutate(across(percent, ~ if_else(is.na(n), "-", as.character(.x)))) %>%
      mutate(across(n, ~ if_else(is.na(.x), "-", scales::comma(.x, accuracy = 1)))) %>%
      mutate(value = as.character(glue("{n} ({percent}%)"))) %>%
      select(arm, !! sym(var), value) %>%
      pivot_wider(
        names_from = arm, 
        values_from = value
      ) %>%
      mutate(variable = var) %>%
      rename("category" = var)
    
    if (is.logical(out$category)) {
      out <- out %>%
        filter(category) %>%
        mutate(across(category, ~ "yes"))
    }
    
    return(out)
    
  }
  
  cat("---- make table 1 ----\n")
  # make table1
  table1 <- bind_rows(lapply(
    variables,
    function(x)
      data %>% summary_var(var = x)
  ))
  
  cat("---- tidy table 1 ----\n")
  # vairables under "History of" heading
  history_of_vars <- c(
    "heart_failure", "other_heart_disease", "dialysis", "diabetes", 
    "chronic_liver_disease", "current_copd", "other_respiratory", "lung_cancer", 
    "haematological_cancer", "cancer_excl_lung_and_haem", "any_immunosuppression", 
    "dementia", "other_neuro_conditions", "ld_inc_ds_and_cp", "psychosis_schiz_bipolar")
  # tidy table1
  table1_tidy <- var_order %>% 
    left_join(var_labels, by = "variable") %>%
    left_join(table1, by = c("category", "variable")) %>%
    mutate(across(category,
                  ~ if_else(variable %in% history_of_vars, variable_label, .x))) %>%
    mutate(across(variable_label,
                  ~ if_else(variable %in% history_of_vars, "History of", .x))) %>%
    mutate(across(variable_label, ~ str_replace(.x, "min_elig_date", as.character(min_elig_date)))) %>%
    rename(Variable = variable_label, Characteristic = category, Unvaccinated = unvax) %>%
    select(-variable) %>%
    select(Variable, Characteristic, everything()) %>%
    mutate(across(c(BNT162b2, ChAdOx, Unvaccinated), 
                  ~ if_else(is.na(.x), "0 (0%)", .x))) 
  
  # age summary
  age_summary <- data %>%
    group_by(arm) %>%
    summarise(
      median = median(age, na.rm=TRUE),
      iqr = IQR(age, na.rm = TRUE),
      .groups = "keep") %>% 
    ungroup() %>%
    transmute(arm, value = glue("{median} ({iqr})")) %>%
    pivot_wider(names_from = "arm", values_from = "value") %>%
    mutate(Variable = "Age", Characteristic = "Median (IQR)") 
  
  age_missing <- data %>% 
    group_by(arm) %>%
    summarise(
      missing = sum(is.na(age)), 
      .groups = "keep") %>%
    ungroup() %>%
    transmute(arm, value = scales::comma(missing, accuracy = 1)) %>%
    pivot_wider(names_from = "arm", values_from = "value") %>%
    mutate(Variable = "Age", Characteristic = "Missing") 
    
  
  table1_tidy_n <- data %>% 
    group_by(arm) %>% 
    count() %>% 
    ungroup() %>% 
    pivot_wider(names_from = arm, values_from = n) %>% 
    mutate(Variable = "", Characteristic = "N") %>%
    mutate(across(c(BNT162b2, ChAdOx, unvax), 
                  ~ scales::comma(.x, accuracy = 1))) %>%
    bind_rows(
      age_summary,
      age_missing
      ) %>%
    rename(Unvaccinated = unvax) %>%
    bind_rows(
      table1_tidy
    )
  
  cat("---- save table1.csv ----\n")
  # save table1_tidy
  readr::write_csv(table1_tidy,
                   here::here("output", "report", "tables", glue("table1_{subgroup_label}_REDACTED.csv")))
  
  cat("---- save table1.html ----\n")
  table1_tidy_n %>%
    gt(
      groupname_col="Variable",
      rowname_col = "Characteristic"
    ) %>%
    tab_header(
      title = glue("Subgroup: {subgroup}"),
      subtitle = "Patient characteristics as of second vaccination period + 2 weeks") %>%
    tab_style(
      style = cell_text(weight="bold"),
      locations = list(
        cells_column_labels(
          columns = everything()
        ),
        cells_row_groups(
          groups = everything()
        ))
    ) %>%
    gtsave(
      filename = glue("table1_{subgroup_label}_REDACTED.html"),
      path = here::here("output", "report", "tables")
    )
  
}






