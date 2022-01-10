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
apply_model_fun <- function(
  subgroup,
  comparison,
  outcome
) {
  splice(
    comment(glue("subgroup = {subgroup}, comparison = {comparison}, outcome = {outcome}")),
    action(
      name = glue("apply_model_cox_{subgroup}_{comparison}_{outcome}"),
      run = "r:latest analysis/comparisons/apply_model_cox.R",
      arguments = c(subgroup, comparison, outcome),
      needs = list(
        "design", 
        "data_comparisons_process", 
        "data_tte_process"),
      highly_sensitive = list(
        modelnumber = glue("output/models/model*_{subgroup}_{comparison}_{outcome}.rds"),
        model_summary_rds = glue("output/models/modelcox_summary_{subgroup}_{comparison}_{outcome}.rds")
      ),
      moderately_sensitive = list(
        incidence_table_all = glue("output/tables/incidence_{subgroup}_{comparison}_{outcome}_all.txt"),
        incidence_table_strata = glue("output/tables/incidence_{subgroup}_{comparison}_{outcome}_strata.txt"),
        model_glance = glue("output/models/modelcox_glance_{subgroup}_{comparison}_{outcome}.csv")
      )
    )
    
  )
}

table_fun <- function(
  subgroup,
  comparison
) {
  
  splice(
    comment(glue("tabulate cox model for all outcomes")),
    action(
      name = glue("tables_model_cox_{subgroup}_{comparison}"),
      run = "r:latest analysis/comparisons/tables_cox.R",
      arguments = c(subgroup, comparison),
      needs = splice("design", 
                     "data_2nd_vax_dates",
                     "data_tte_process",
                     lapply(
                       outcomes, 
                       function(x) 
                         glue("apply_model_cox_{subgroup}_{comparison}_{x}"))),
      moderately_sensitive = list(
        table_glance = glue("output/tables/modelcox_glance_{subgroup}_{comparison}.txt"),
        table_coefficients = glue("output/tables/modelcox_coefficients_{subgroup}_{comparison}.txt"))
    )
  )
  
  
}

plot_fun <- function(
  plot
) {
  
  if (plot %in% c("BNT162b2", "ChAdOx")) {
    comparisons <- plot
  } else if (plot %in% "BNT162b2andChAdOx") {
    comparisons <- c("BNT162b2", "ChAdOx")
  } else if (plot %in% "BNT162b2vsChAdOx") {
    comparisons <- "both"
  }
  
  if (str_detect(plot, "ChAdOx")) {
    subgroups <- "03-10"
  }
  
  splice(
    comment(glue("plot {plot}")),
    action(
      name = glue("plot_model_cox_{plot}"),
      run = "r:latest analysis/comparisons/plot_cox.R",
      arguments = plot,
      needs = splice("design",
                     "data_2nd_vax_dates",
                     as.list(unlist(lapply(
                       subgroups,
                       function(x)
                         unlist(lapply(
                           comparisons,
                           function(y)
                             unlist(lapply(
                               outcomes,
                               function(z)
                                 glue("apply_model_cox_{x}_{y}_{z}")
                             ), recursive = FALSE)
                         ), recursive = FALSE)
                     ), recursive = FALSE))),
      moderately_sensitive = list(
        plot = glue("output/images/plot_res_{plot}.png"))
    )
  )
  
  
}
    
# specify project ----

## defaults ----
defaults_list <- list(
  version = "3.0",
  expectations= list(population_size=100000L)
)



outcomes <- c("postest", "covidadmitted", "coviddeath", "death")
subgroups <- c("02", "03-10", "11-12")  
comparisons <- c("BNT162b2", "ChAdOx", "both")
plots <- c("BNT162b2", "ChAdOx", "BNT162b2andChAdOx", "BNT162b2vsChAdOx")

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
      data_vax_dates = "output/data/data_*_vax_dates.rds",
      data_long_dates = "output/data/data_long_*_dates.rds",
      data_covid_any = "output/data/data_covid_any.rds"
    ),
    moderately_sensitive = list(
      data_properties = "output/tables/data_processed_tabulate.txt"
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
      second_vax_period_dates_txt = "output/tables/second_vax_period_dates.txt"
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
  
  comment("####################################",
          "comparisons", 
          "####################################"),
  
  comment("process comparisons data"),
  action(
    name = "data_comparisons_process",
    run = "r:latest analysis/comparisons/data_comparisons_process.R",
    needs = list(
      "design", 
      "data_input_process", 
      "data_2nd_vax_dates", 
      "data_eligible_cd"),
    highly_sensitive = list(
      data_comparisons = glue("output/data/data_comparisons_*.rds")
    )
  ),
  
  comment(glue("process outcomes data")),
  action(
    name = "data_outcomes_process",
    run = "r:latest analysis/comparisons/data_outcomes_process.R",
    needs = list(
      "design",
      "data_input_process",
      "data_comparisons_process"),
    highly_sensitive = list(
      data_outcomes = "output/data/data_outcomes_*.rds"
    )
  ),
  
  comment(glue("process tte data")),
  action(
    name = "data_tte_process",
    run = "r:latest analysis/comparisons/data_tte_process.R",
    needs = list(
      "design",
      "data_comparisons_process", 
      "data_outcomes_process"),
    highly_sensitive = list(
      data_tte_brand_outcome = "output/data/data_tte_*.rds"
    )
  ),
  
  comment("apply models"),
  splice(
    # over subgroups
    unlist(lapply(
      subgroups, # subgroups
      function(x) {
        
        if (x %in% c("02", "11-12")) {
          comparisons <- "BNT162b2"
        } 
        
        # over comparisons
        unlist(lapply(
          comparisons,
          function(y)
            unlist(lapply(
              outcomes,
              function(z)
              apply_model_fun(
                subgroup = x,
                comparison = y, 
                outcome = z)  
            ),
            recursive = FALSE)
        ),
        recursive = FALSE)
      }
      
    ), recursive = FALSE)
    #
  ),
  
  comment("generate tables"),
  splice(
    # over subgroups
    unlist(lapply(
      subgroups, # subgroups
      function(x) {
        
        if (x %in% c("02", "11-12")) {
          comparisons <- "BNT162b2"
        } 
        
        # over comparisons
        unlist(lapply(
          comparisons,
          function(y)
                table_fun(
                  subgroup = x,
                  comparison = y
                  )  
        ),
        recursive = FALSE)
      }
      
    ), recursive = FALSE)
    #
  ),
  
  comment("generate plots"),
  splice(
    unlist(lapply(plots,
                  function(p)
                    plot_fun(plot = p)
                  ), recursive = FALSE))
  
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
