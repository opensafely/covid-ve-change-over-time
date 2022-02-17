################################################################################

# This script:
# applies the cox model for a given comparison and outcome


################################################################################
library(tidyverse)
library(glue)
library(survival)
library(multcomp)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "BNT162b2"
  subgroup_label <- 1
  outcome <- "anytest"
  
} else{
  comparison <- args[[1]]
  subgroup_label <- as.integer(args[[2]])
  outcome <- args[[3]]
}

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))
K <- study_parameters$max_comparisons

# specify filename_suffix for saving models
filename_suffix <- glue("{comparison}_{subgroup_label}_{outcome}")

# read adjusted model
model2 <- readr::read_rds(
  here::here("output", "models_cox", "data", glue("model2_{filename_suffix}.rds")))

################################################################################


p <- length(model2$coefficients)

names_ref_1 <- sapply(2:K, function(x) glue("k{x} - k1"))
mat_ref_1 <- matrix(0, nrow=K, ncol=p)
mat_ref_1[,1] <- -1
for (i in 2:K) {
  mat_ref_1[i,i] <- 1
}


names_ref_prev <- sapply(2:K, function(x) glue("k{x} - k{x-1}"))
mat_ref_prev <- matrix(0, nrow=K, ncol=p)
mat_ref_prev[,1] <- -1
for (i in 2:K) {
  mat_ref1[i,i] <- 1
}




lincom_ref1 <- matrix(0, nrow = )


summary(model1_BNT162b2_1_anytest)

f <- model_input[[2]]$unadjusted


mod2 <- coxph(Surv(tstart, tstop, status, type = "counting") ~ relevel(k, ref = "1") + strata(strata_var), data = data_cox)


lincom <- rbind("k2 - k1" = c(-1, 1, 0, 0, 0, 0),     
                "k3 - k1" = c(-1, 0, 1, 0, 0, 0))
ref1 <- glht(
  model1_BNT162b2_1_anytest,
  linfct = lincom
)
sref1 <- summary(ref1)
