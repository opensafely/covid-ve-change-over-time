################################################################################

# This script:
# applies the cox model for a given comparison and outcome


################################################################################
library(tidyverse)
library(glue)
library(survival)
library(broom)
library(broom.helpers)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "BNT162b2"
  subgroup_label <- 4
  outcome <- "covidadmitted"
  
} else{
  comparison <- args[[1]]
  subgroup_label <- as.integer(args[[2]])
  outcome <- args[[3]]
}

################################################################################
# create directories
fs::dir_create(here::here("output", "models_cox", "data"))
fs::dir_create(here::here("output", "models_cox", "tables"))

################################################################################
# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroup <- subgroups[subgroup_label]

# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))
K <- study_parameters$max_comparisons

################################################################################
for (kk in 1:K) {
  
  # read data
  model_input <-  readr::read_rds(
    here::here("output", "preflight", "data", glue("model_input_{comparison}_{subgroup_label}_{outcome}_{kk}.rds"))
  )
  
  # specify filename_suffix for saving models
  filename_suffix <- glue("{comparison}_{subgroup_label}_{outcome}_{kk}")
  
  
  if (is.null(model_input)) {
    # model input null if not enough events, see preflight.R for threshold
    # save empty tibbles to avoid errors
    readr::write_rds(
      tibble(),
      here::here("output", "models_cox", "data", glue("modelcox_glance_{filename_suffix}.rds"))) 
    
    readr::write_rds(
      tibble(),
      here::here("output", "models_cox", "data", glue("modelcox_tidy_{filename_suffix}.rds"))) 
    
    readr::write_rds(
      tibble(),
      here::here("output", "models_cox", "data", glue("modelx_{comparison}_{filename_suffix}.rds")), 
      compress="gz")
    
  } else {
    
    model_varlist <- readr::read_rds(
      here::here("output", "lib", "model_varlist.rds")
    )
    
    data_cox <- model_input$data
    formulas <- model_input$formulas
    
    ################################################################################
    # define formulas
    
    model_names = c(
      "1" = "Unadjusted",
      "2" = "Adjusting for demographics + clinical"
    )
    
    formula_cox_1 <- formulas$unadjusted
    formula_cox_2 <- formulas$unadjusted %>% 
      update(formulas$demographic) %>% 
      update(formulas$clinical)
    
    ################################################################################
    # apply and save models (coxph within functions was causing dramas with broom)
    
    # set max iterations
    opt_control <- coxph.control(iter.max = 30)
    
    coxmods <- list()
    cat(glue("...... fitting unadjusted model ......"), "\n")
    coxmods[[1]] <- coxph(
      formula = formula_cox_1,
      data = data_cox,
      # robust = TRUE,
      # id = patient_id,
      cluster = patient_id,
      na.action = "na.fail",
      control = opt_control)
    
    readr::write_rds(
      coxmods[[1]],
      here::here("output", "models_cox", "data", glue("model1_{filename_suffix}.rds")),
      compress = "gz"
    )
    
    cat(glue("...... fitting adjusted model ......"), "\n")
    coxmods[[2]] <- coxph(
      formula = formula_cox_2,
      data = data_cox,
      # robust = TRUE,
      # id = patient_id,
      cluster = patient_id,
      na.action = "na.fail",
      control = opt_control)
    
    readr::write_rds(
      coxmods[[2]],
      here::here("output", "models_cox", "data", glue("model2_{filename_suffix}.rds")),
      compress = "gz"
    )
    
    ################################################################################
    # process output
    glance <- list()
    tidy <- list()
    for (i in seq_along(coxmods)) {
      
      coxmod <- coxmods[[i]]
      
      glance[[i]] <-
        broom::glance(coxmod) %>%
        add_column(
          model = i,
          convergence = coxmod$info[["convergence"]],
          ram = format(object.size(coxmod), units="GB", standard="SI", digits=3L),
          .before = 1
        )
      
      tidy[[i]] <-
        broom.helpers::tidy_plus_plus(
          coxmod,
          exponentiate = FALSE
        ) %>%
        tidy_add_n(
          coxmod
        ) %>%
        add_column(
          model = i,
          .before=1
        )
      
    }
    
    model_glance <- bind_rows(glance)
    readr::write_rds(
      model_glance,
      here::here("output", "models_cox", "data", glue("modelcox_glance_{filename_suffix}.rds"))) 
    
    model_tidy <- bind_rows(tidy)
    readr::write_rds(
      model_tidy,
      here::here("output", "models_cox", "data", glue("modelcox_tidy_{filename_suffix}.rds"))) 
  }
  
}







