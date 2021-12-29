################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobs <- FALSE
  group <- "02"
  outcome <- "postest"
  
} else{
  removeobs <- TRUE
  group <- args[[1]]
  outcome <- args[[2]]
}

################################################################################

fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "tables"))