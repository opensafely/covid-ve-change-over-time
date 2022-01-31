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
  subgroup_label,
  comparison,
  outcome
) {
  
  subgroup <- subgroups[subgroup_label]
  
  splice(
    comment(glue("{comparison}; {subgroup}; {outcome}")),
    comment("preflight checks"),
    action(
      name = glue("preflight_{comparison}_{subgroup_label}_{outcome}"),
      run = "r:latest analysis/comparisons/preflight.R",
      arguments = c(comparison, subgroup_label, outcome),
      needs = list(
        "design", 
        "data_comparisons_process", 
        glue("data_tte_process_{comparison}")
        ),
      highly_sensitive = list(
        model_input = glue("output/preflight/data/model_input_{comparison}_{subgroup_label}_{outcome}.rds")
      ),
      moderately_sensitive = list(
        eventcheck_table = glue("output/preflight/tables/eventcheck_{comparison}_{subgroup_label}_{outcome}_*.html"),
        preflight_report = glue("output/preflight/tables/preflight_report_{comparison}_{subgroup_label}_{outcome}.txt")
      )
    ),
    comment("apply cox model"),
    action(
      name = glue("apply_model_cox_{comparison}_{subgroup_label}_{outcome}"),
      run = "r:latest analysis/comparisons/apply_model_cox.R",
      arguments = c(comparison, subgroup_label, outcome),
      needs = list(
        "design", 
        glue("preflight_{comparison}_{subgroup_label}_{outcome}")),
      highly_sensitive = list(
        modelnumber = glue("output/models_cox/data/model*_{comparison}_{subgroup_label}_{outcome}.rds"),
        model_summary = glue("output/models_cox/data/modelcox_summary_{comparison}_{subgroup_label}_{outcome}.rds"),
        model_glance = glue("output/models_cox/data/modelcox_glance_{comparison}_{subgroup_label}_{outcome}.rds")
      )
    )
    
  )
}

table_fun <- function(
  comparison,
  subgroup_label
) {
  
  splice(
    comment(glue("tabulate cox model for all outcomes")),
    action(
      name = glue("tables_model_cox_{comparison}_{subgroup_label}"),
      run = "r:latest analysis/comparisons/tables_cox.R",
      arguments = c(comparison, subgroup_label),
      needs = splice("design", 
                     "data_2nd_vax_dates",
                     glue("data_tte_process_{comparison}"),
                     lapply(
                       outcomes_model, 
                       function(x) 
                         glue("apply_model_cox_{comparison}_{subgroup_label}_{x}"))),
      moderately_sensitive = list(
        table_glance = glue("output/models_cox/tables/modelcox_glance_{comparison}_{subgroup_label}.txt"),
        table_coefficients = glue("output/models_cox/tables/modelcox_coefficients_{comparison}_{subgroup_label}.txt"))
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
    subgroup_labels <- subgroup_labels[-which(subgroups == "18-39")]
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
                       comparisons,
                       function(x)
                         unlist(lapply(
                           subgroup_labels,
                           function(y)
                             unlist(lapply(
                               outcomes_model,
                               function(z)
                                 glue("apply_model_cox_{x}_{y}_{z}")
                             ), recursive = FALSE)
                         ), recursive = FALSE)
                     ), recursive = FALSE))),
      moderately_sensitive = list(
        plot = glue("output/models_cox/images/plot_res_{plot}_*.png"))
    )
    
  )
  
  
}
    
# specify project ----

## defaults ----
defaults_list <- list(
  version = "3.0",
  expectations= list(population_size=100000L)
)



outcomes <- readr::read_rds(here::here("output", "lib", "outcomes.rds"))
subgroups <- c(readr::read_rds(here::here("output", "lib", "subgroups.rds")), "all")
subgroup_labels <- seq_along(subgroups)
comparisons <- c("BNT162b2", "ChAdOx", "both")
plots <- c("BNT162b2", "ChAdOx", "BNT162b2andChAdOx", "BNT162b2vsChAdOx")
outcomes_model <- outcomes#[-which(outcomes=="noncoviddeath")]

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
      outcomes = "output/lib/outcomes.rds",
      subgroups = "output/lib/subgroups.rds"
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
      data_all = "output/data/data_*.rds"
      # data_processed = "output/data/data_processed.rds",
      # data_vax_dates = "output/data/data_*_vax_dates.rds",
      # data_long_dates = "output/data/data_long_*_dates.rds",
      # data_covid_any = "output/data/data_covid_any.rds"
    ),
    moderately_sensitive = list(
      data_properties = "output/tables/data_*_tabulate.txt"
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
      second_vax_period_dates_rds = "output/second_vax_period/data/second_vax_period_dates.rds"
    ),
    moderately_sensitive = list(
      second_vax_period_dates_txt = "output/second_vax_period/tables/second_vax_period_dates.txt"
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
  
  comment("apply eligiblity criteria from boxes c, d and e"),
  action(
    name = "data_eligible_cde",
    run = "r:latest analysis/second_vax_period/data_eligible_cde.R",
    needs = list("design", "data_input_process", "data_eligible_ab", "data_2nd_vax_dates"),
    highly_sensitive = list(
      data_eligible_e_vax = "output/data/data_eligible_e_vax.rds",
      data_eligible_e_unvax = "output/data/data_eligible_e_unvax.rds",
      data_eligible_e = "output/data/data_eligible_e.csv"
    )
  ),
  
  comment("####################################", 
          "study definition tests",
          "####################################"),
  # comment("generate dummy data for study_definition_tests"),
  # action(
  #   name = "dummy_data",
  #   run = "r:latest analysis/dummy_data.R",
  #   needs = list("design"),
  #   moderately_sensitive = list(
  #     dummy_data = "analysis/dummy_data.feather"
  #   )
  # ),
  
  comment("study definition tests"),
  action(
    name = "generate_covid_tests_data",
    run = "cohortextractor:latest generate_cohort --study-definition study_definition_tests --output-format feather",
    # dummy_data_file = "analysis/dummy_data.feather",
    needs = list("design", "data_eligible_cde"),
    highly_sensitive = list(
      cohort = "output/input_tests.feather"
    )
  ),
  
  comment("check the tests data as expected"),
  action(
    name = "check_tests",
    run = "r:latest analysis/tests/check_tests.R",
    needs = list("design", "generate_covid_tests_data"),
    moderately_sensitive = list(
      covariate_distribution = "output/tests/images/covariate_distribution.png",
      data_tests_tabulate = "output/tests/tables/data_tests_tabulate.txt"
    )
  ),
  
  comment("####################################",
          "subsequent vaccination", 
          "####################################"),
  
  comment("plot cumulative incidence of subsequent vaccination"),
  action(
    name = "plot_cumulative_incidence",
    run = "r:latest analysis/subsequent_vax/plot_cumulative_incidence.R",
    needs = list("design", "data_input_process", "data_eligible_cde"),
    moderately_sensitive = list(
      ci_vax = "output/subsequent_vax/images/ci_vax_*.png",
      survtable = "output/subsequent_vax/tables/survtable_*.txt"
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
      "data_eligible_cde"),
    highly_sensitive = list(
      data_comparisons = glue("output/comparisons/data/data_comparisons_*.rds")
    )
  ),
  
  comment(glue("process tte data")),
  splice(unlist(lapply(
    comparisons,
    function(x)
      action(
        name = glue("data_tte_process_{x}"),
        run = "r:latest analysis/comparisons/data_tte_process.R",
        arguments = x,
        needs = list(
          "design",
          "data_input_process",
          "data_comparisons_process"),
        highly_sensitive = list(
          data_tte_brand_outcome = glue("output/tte/data/data_tte_{x}*.rds"),
          event_counts = glue("output/tte/tables/event_counts_{x}.rds")
        )
      )
  ), recursive = FALSE)),
  
  comment(glue("process event counts tables")),
  action(
    name = glue("process_event_count_tables"),
    run = "r:latest analysis/comparisons/process_event_count_tables.R",
    needs = list(
      "design",
      "data_tte_process_BNT162b2",
      "data_tte_process_ChAdOx"),
    moderately_sensitive = list(
      tidy_tables_events = glue("output/tte/tables/tidy_events*.txt")
    )
  ),
  
  comment("apply models"),
  splice(
    # over subgroups
    unlist(lapply(
      comparisons,

      function(x) {
        if (!(x %in% "BNT162b2")) {
          subgroup_labels <- subgroup_labels[-which(subgroups == "18-39")]
        }
        unlist(lapply(
          subgroup_labels,

          function(y)
            splice(
            unlist(lapply(
              outcomes_model,

              function(z)
              apply_model_fun(
                comparison = x,
                subgroup_label = y,
                outcome = z)
            ),
            recursive = FALSE),
            table_fun(
              comparison = x,
              subgroup_label = y
            )
            )
          
        ),
        recursive = FALSE)
      }
    ), recursive = FALSE)
  ),

  comment("generate plots"),
  splice(
    unlist(lapply(plots,
                  function(p)
                    plot_fun(plot = p)
                  ), recursive = FALSE))#,
  # 
  # comment("combine all incidence tables"),
  # action(
  #   name = glue("combine_incidence_tables"),
  #   run = "r:latest analysis/comparisons/combine_incidence_tables.R",
  #   needs = splice("design", 
  #                  # "data_2nd_vax_dates",
  #                  # "data_tte_process",
  #                  as.list(unlist(lapply(
  #                    subgroups, # subgroups
  #                    function(x) {
  #                      
  #                      if (x %in% c("02", "11-12")) {
  #                        comparisons <- "BNT162b2"
  #                      } 
  #                      
  #                      # over comparisons
  #                      unlist(lapply(
  #                        comparisons,
  #                        function(y)
  #                          unlist(lapply(
  #                            outcomes,
  #                            function(z)
  #                              glue("apply_model_cox_{x}_{y}_{z}")
  #                          ),
  #                          recursive = FALSE)
  #                      ),
  #                      recursive = FALSE)
  #                    }
  #                    
  #                  ), recursive = FALSE))),
  #   moderately_sensitive = list(
  #     incidence_all = glue("output/tables/incidence_all.txt"))
  # ),
  # 
  # comment("plot cumulative incidence of 1st or 3rd vacciantion"),
  # splice(unlist(lapply(
  #   c("BNT162b2", "ChAdOx"),
  #   function(x)
  #     action(
  #       name = glue("plot_cumulative_incidence_vax_{x}"),
  #       run = "r:latest analysis/comparisons/plot_cumulative_incidence_vax.R",
  #       arguments = x,
  #       needs = list(
  #         "design",
  #         "data_input_process",
  #         "data_comparisons_process"),
  #       moderately_sensitive = list(
  #         plot_cumulative_incidence_vax = glue("output/images/cumulative_incidence_{x}.png"),
  #         survtable = glue("output/tables/survtable_{x}.txt")
  #       )
  #     )
  # ), recursive = FALSE))
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
