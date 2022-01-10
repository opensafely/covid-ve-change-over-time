library(tidyverse)
library(glue)
library(kableExtra)

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))

process_table <- function(
  subgroup,
  comparison,
  outcome
) {
    
    table <- read_delim(
      here::here("output", "tables", glue("incidence_{subgroup}_{comparison}_{outcome}_all.txt")), 
      delim = "|", escape_double = FALSE, col_names = FALSE, 
      trim_ws = TRUE, skip = 4, progress = FALSE)
    
    table_col_names <- unname(unlist(table[1,-c(1,5)]))
    
    table[-c(1:3), -c(1,5)] %>%
      kable("pipe",
            col.names  = table_col_names,
            caption = str_c(glue("---- subgroup: {subgroup}; comparison: {comparison}; outcome: {outcome} ----\nevents / person-years"),
                            " (n_{k}, percent of n_{k-1} for vax and n_{k-2} for unvax)"))
    
}

comparisons <- c("BNT162b2", "ChAdOx")


file.remove(
  here::here("output", "tables", glue("incidence_all.txt"))
)
for (y in comparisons) {
  
  if (y == "BNT162b2") subgroups <- c("02", "03-10", "11-12") else subgroups <- "03-10"
  
  for (x in subgroups) {
    for (z in outcomes) {
      
      capture.output(
        process_table(subgroup=x,comparison=y,outcome=z),
        file = here::here("output", "tables", glue("incidence_all.txt")),
        append=TRUE
      )
      
    }
  }
}
