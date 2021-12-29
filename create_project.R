library('tidyverse')
library('yaml')
library('here')
library('glue')

# create action functions ----

## generic action function ----
action <- function(
  name,
  run,
  dummy_data_file=NULL,
  arguments=NULL,
  needs=NULL,
  highly_sensitive=NULL,
  moderately_sensitive=NULL
){
  
  outputs <- list(
    highly_sensitive = highly_sensitive,
    moderately_sensitive = moderately_sensitive
  )
  outputs[sapply(outputs, is.null)] <- NULL
  
  action <- list(
    run = paste(c(run, arguments), collapse=" "),
    dummy_data_file = dummy_data_file,
    needs = needs,
    outputs = outputs
  )
  action[sapply(action, is.null)] <- NULL
  
  action_list <- list(name = action)
  names(action_list) <- name
  
  action_list
}


## create comment function ----
comment <- function(...){
  list_comments <- list(...)
  comments <- map(list_comments, ~paste0("## ", ., " ##"))
  comments
}


## create function to convert comment "actions" in a yaml string into proper comments
convert_comment_actions <-function(yaml.txt){
  yaml.txt %>%
    str_replace_all("\\\n(\\s*)\\'\\'\\:(\\s*)\\'", "\n\\1")  %>%
    #str_replace_all("\\\n(\\s*)\\'", "\n\\1") %>%
    str_replace_all("([^\\'])\\\n(\\s*)\\#\\#", "\\1\n\n\\2\\#\\#") %>%
    str_replace_all("\\#\\#\\'\\\n", "\n")
}

## actions that extract and process data ----

actions_comparisons <- function(
  jcvi_group, outcomes, # arguments
  check_combine_outcomes=TRUE,
  apply_models=TRUE
  ) {
  
  actions <- splice(
    
    comment("####################################",
            "comparisons", 
            "####################################"),
    comment(glue("process data, apply model and generate report for JCVI group {jcvi_group}")),
    comment("process covariates data"),
    action(
      name = glue("data_comparisons_process_{jcvi_group}"),
      run = "r:latest analysis/comparisons/data_comparisons_process.R",
      arguments = c(jcvi_group),
      needs = list("design", "data_input_process", "data_long_process", "data_2nd_vax_dates", "data_eligible_cd"),
      highly_sensitive = list(
        data_comparisons = glue("output/jcvi_group_{jcvi_group}/data/data_comparisons.rds")
      )
    ),
    
    comment(glue("process outcomes data")),
    action(
      name = glue("data_outcomes_process_{jcvi_group}"),
      run = "r:latest analysis/comparisons/data_outcomes_process.R",
      arguments = c(jcvi_group),
      needs = list("design", "data_input_process", "data_long_process", glue("data_comparisons_process_{jcvi_group}")),
      highly_sensitive = list(
        data_outcomes = glue("output/jcvi_group_{jcvi_group}/data/data_outcomes.rds")
      )
    )
  )
  
  if (check_combine_outcomes) {
    actions <- splice(
      actions,
      
      comment(glue("check gap between outcomes for combining")),
      action(
        name = glue("check_combine_outcomes_{jcvi_group}"),
        run = "r:latest analysis/comparisons/check_combine_outcomes.R",
        arguments = c(jcvi_group),
        needs = list(glue("data_outcomes_process_{jcvi_group}")),
        highly_sensitive = list(
          data_check_combine_outcomes = glue("output/jcvi_group_{jcvi_group}/data/check_combine_outcomes.rds")
        ),
        moderately_sensitive = list(
          plot_check_combine_outcomes = glue("output/jcvi_group_{jcvi_group}/images/check_combine_outcomes.png"),
          table_check_combine_outcomes = glue("output/jcvi_group_{jcvi_group}/tables/check_combine_outcomes.csv")
        )
      )
    )
  }
  
  # define function for apply models actions
  if (apply_models) {
    
    apply_models_fun <- function(outcome) {
      splice(
        
        comment(glue("outcome = {outcome}")),
        comment(glue("process tte data for {outcome}")),
        action(
          name = glue("data_tte_process_{jcvi_group}_{outcome}"),
          run = "r:latest analysis/comparisons/data_tte_process.R",
          arguments = c(jcvi_group, outcome),
          needs = list(glue("data_comparisons_process_{jcvi_group}"), glue("data_outcomes_process_{jcvi_group}")),
          highly_sensitive = list(
            data_tte_brand_outcome = glue("output/jcvi_group_{jcvi_group}/data/data_tte_*_{outcome}.rds")
          ),
        ),
        
        comment(glue("apply cox model for {outcome}")),
        action(
          name = glue("apply_model_cox_{jcvi_group}_{outcome}"),
          run = "r:latest analysis/comparisons/apply_model_cox.R",
          arguments = c(jcvi_group, outcome),
          needs = list("design", glue("data_comparisons_process_{jcvi_group}"), glue("data_tte_process_{jcvi_group}_{outcome}")),
          highly_sensitive = list(
            modelnumber = glue("output/jcvi_group_{jcvi_group}/models/*_{outcome}_model*.rds"),
            model_tidy_rds = glue("output/jcvi_group_{jcvi_group}/models/*_{outcome}_modelcox_tidy.rds"),
            model_summary_rds = glue("output/jcvi_group_{jcvi_group}/models/*_{outcome}_modelcox_summary.rds")
          ),
          moderately_sensitive = list(
            model_glance = glue("output/jcvi_group_{jcvi_group}/models/*_{outcome}_modelcox_glance.csv"),
            model_tidy_csv = glue("output/jcvi_group_{jcvi_group}/models/*_{outcome}_modelcox_tidy.csv")
          )
        ),
        
        comment(glue("plot cox model for all outcomes")),
        action(
          name = glue("plot_model_cox_{jcvi_group}"),
          run = "r:latest analysis/comparisons/plot_cox.R",
          arguments = c(jcvi_group),
          needs = splice("design", "data_2nd_vax_dates",
                         lapply(
                           outcomes, 
                           function(x) glue("apply_model_cox_{jcvi_group}_{x}"))),
          moderately_sensitive = list(
            plot = glue("output/jcvi_group_{group}/images/plot_res_*.png"))
          ),
        
        comment(glue("tabulate cox model for all outcomes")),
        action(
          name = glue("tables_model_cox_{jcvi_group}"),
          run = "r:latest analysis/comparisons/tables_cox.R",
          arguments = c(jcvi_group),
          needs = splice("design", "data_2nd_vax_dates",
                         lapply(
                           outcomes, 
                           function(x) glue("apply_model_cox_{jcvi_group}_{x}"))),
          moderately_sensitive = list(
            table_glance = glue("output/jcvi_group_{group}/tables/*_modelcox_glance.txt"),
            table_coefficients = glue("output/jcvi_group_{group}/tables/*_modelcox_coefficients.txt")
            )
        )
        
      )
    }
    
    
    actions <- splice(
      
      actions,
      # apply models actions for all outcomes
      unlist(lapply(
        outcomes,
        apply_models_fun
      ),
      recursive = FALSE
      )
    )
  }
  
  return(actions)
  
}

# specify project ----

## defaults ----
defaults_list <- list(
  version = "3.0",
  expectations= list(population_size=100000L)
)

## actions ----
actions_list <- splice(
  
  comment("####################################",
          "preliminaries",
          "####################################"),
  action(
    name = "design",
    run = "r:latest analysis/design.R",
    moderately_sensitive = list(
      study_dates_json = "output/lib/study_parameters.json",
      study_dates_rds = "output/lib/study_parameters.rds",
      jcvi_groups = "output/lib/jcvi_groups.csv",
      elig_dates = "output/lib/elig_dates.csv",
      regions = "output/lib/regions.csv",
      model_varlist = "output/lib/model_varlist.rds",
      outcomes = "output/lib/outcomes.rds"
    )
  ),
  
  comment("####################################", 
          "study definition",
          "####################################"),
  comment("generate dummy data for study_definition"),
  action(
    name = "dummy_data",
    run = "r:latest analysis/dummy_data.R",
    needs = list("design"),
    moderately_sensitive = list(
      dummy_data = "analysis/dummy_data.feather"
    )
  ),
  
  comment("study definition"),
  action(
    name = "generate_study_population",
    run = "cohortextractor:latest generate_cohort --study-definition study_definition --output-format feather",
    dummy_data_file = "analysis/dummy_data.feather",
    needs = list("design", "dummy_data"),
    highly_sensitive = list(
      cohort = "output/input.feather"
    )
  ),
  
  comment("####################################", 
          "preprocessing",
          "####################################"),
  
  comment("process data from study_definition"),
  action(
    name = "data_input_process",
    run = "r:latest analysis/preprocess/data_input_process.R",
    needs = list("design", "dummy_data", "generate_study_population"),
    highly_sensitive = list(
      data_covs = "output/data/data_covs.rds",
      data_vax_dates = "output/data/data_*_vax_dates.rds"
    ),
    moderately_sensitive = list(
      data_properties = "output/tables/data_processed_tabulate.txt"
    )
  ),
  
  comment("process recurring variables as long data"),
  action(
    name = "data_long_process",
    run = "r:latest analysis/preprocess/data_long_process.R",
    needs = list("design", "data_input_process"),
    highly_sensitive = list(
      data_long_dates = "output/data/data_long_*_dates.rds"
    )
  ),
  
  comment("apply eligiblity criteria from boxes a and b"),
  action(
    name = "data_eligible_ab",
    run = "r:latest analysis/preprocess/data_eligible_ab.R",
    needs = list("design", "data_input_process"),
    highly_sensitive = list(
      data_eligible_a = "output/data/data_eligible_a.rds",
      data_eligible_b = "output/data/data_eligible_b.rds"
    ),
    moderately_sensitive = list(
      eligibility_count = "output/lib/eligibility_count_ab.csv",
      group_age_ranges = "output/lib/group_age_ranges.csv"
    )
  ),
  
  comment("####################################", 
          "second_vax_period",
          "####################################"),
  comment("identify second vaccination time periods"),
  comment("create dataset for identifying second vaccination time periods"),
  action(
    name = "data_2nd_vax_dates",
    run = "r:latest analysis/second_vax_period/data_2nd_vax_dates.R",
    needs = list("design", "data_input_process", "data_eligible_ab"),
    highly_sensitive = list(
      data_vax_plot = "output/second_vax_period/data/data_vax_plot.rds",
      second_vax_period_dates_rds = "output/lib/second_vax_period_dates.rds"
    ),
    moderately_sensitive = list(
      second_vax_period_dates_csv = "output/lib/second_vax_period_dates.csv"
    )
  ),
  
  comment("plot second vaccination time periods"),
  action(
    name = "plot_2nd_vax_dates",
    run = "r:latest analysis/second_vax_period/plot_2nd_vax_dates.R",
    needs = list("design", "data_eligible_ab", "data_2nd_vax_dates"),
    moderately_sensitive = list(
      plots_by_region = "output/second_vax_period/images/plot_by_region_*.png"
    )
  ),
  
  comment("apply eligiblity criteria from boxes c and d"),
  action(
    name = "data_eligible_cd",
    run = "r:latest analysis/second_vax_period/data_eligible_cd.R",
    needs = list("design", "data_input_process", "data_eligible_ab", "data_2nd_vax_dates"),
    highly_sensitive = list(
      data_eligible_c = "output/data/data_eligible_c.rds",
      data_eligible_d = "output/data/data_eligible_d.rds"
    )
  ),
  
  actions_comparisons(jcvi_group = "02", outcomes = "postest")
  
)

## combine everything ----
project_list <- splice(
  defaults_list,
  list(actions = actions_list)
)

## convert list to yaml, reformat comments and whitespace,and output ----
as.yaml(project_list, indent=2) %>%
  # convert comment actions to comments
  convert_comment_actions() %>%
  # add one blank line before level 1 and level 2 keys
  str_replace_all("\\\n(\\w)", "\n\n\\1") %>%
  str_replace_all("\\\n\\s\\s(\\w)", "\n\n  \\1") %>%
  writeLines(here("project.yaml"))


## grab all action names and send to a txt file

names(actions_list) %>% tibble(action=.) %>%
  mutate(
    model = str_detect(action, "model"),
    model_number = cumsum(model)
  ) %>%
  group_by(model_number) %>%
  summarise(
    sets = paste(action, collapse=" ")
  ) %>% pull(sets) %>%
  paste(collapse="\n") %>%
  writeLines(here("actions.txt"))
