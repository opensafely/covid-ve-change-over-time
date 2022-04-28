library(tidyverse)
library(glue)
library(kableExtra)

################################################################################
if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)


################################################################################
# prepare the table 1 data
table_out0 <- bind_rows(
  lapply(
    1:4,
    function(x) 
      readr::read_csv(here::here(release_folder, glue("table1_{x}_REDACTED.csv")),
                      show_col_types = FALSE) %>%
      mutate(subgroup = x)
  ) 
) %>%
  pivot_wider(
    names_from = subgroup,
    values_from = c(BNT162b2, ChAdOx1, Unvaccinated),
    names_glue = "{subgroup}_{.value}"
  ) %>%
  select(-`4_ChAdOx1`) %>%
  select(Variable, Characteristic, starts_with(as.character(subgroup_labels))) %>%
  mutate(across(starts_with(as.character(subgroup_labels)),
                ~ if_else(
                  is.na(.x),
                  "-",
                  .x)))

################################################################################
# prepare the column names
# column names for table
table_names <- names(table_out0)
# remove subgroup label
col_names <- str_remove(table_names, "\\d_")
names(col_names) <- table_names
# number of columns for each subgroup
cols_subtype <- sapply(
  as.character(subgroup_labels), 
  function(x) sum(str_detect(table_names, glue("{x}_")))
)
names(cols_subtype) <- subgroups
# reorder
cols_subtype <- cols_subtype[subgroup_labels]

################################################################################
# create and save the version for the manuscript
variable_order_manuscript <- c(
  NA_character_, 
  "Age", 
  "Sex", 
  "IMD", 
  "Ethnicity",
  "BMI", 
  "Morbidity count",
  "Number of SARS-CoV-2 tests between 2020-05-18 and 2020-12-08", 
  "Flu vaccine in previous 5 years")

table1_data_manuscript <- tibble(Variable = variable_order_manuscript) %>%
  left_join(table_out0, by = "Variable") %>%
  mutate(across(Variable, ~str_remove(.x, " between 2020-05-18 and 2020-12-08"))) %>%
  mutate(across(Variable, ~str_remove(.x, " criteria met"))) 

table1_manuscript <- list(
  col_names = col_names,
  cols_subtype = cols_subtype,
  data = table1_data_manuscript
)

readr::write_rds(
  table1_manuscript,
  here::here(release_folder, "table1_manuscript.rds"))

################################################################################
# create and save the version for the supplement

# add footnote marker to not clinically vulnerable subgroups for supplement
names(cols_subtype)[3:4] <- str_c(
  names(cols_subtype)[3:4],
  footnote_marker_alphabet(1, format = "latex", double_escape = TRUE)
)

# create the version for the appendix
variable_order_supplement <- c(
  NA_character_,
  "Region", 
  "JCVI group", 
  "Evidence of", 
  "Pregnancy")

table1_data_supplement <- tibble(Variable = variable_order_supplement) %>%
  left_join(table_out0, by = "Variable") %>%
  mutate(tmp = row_number()) %>%
  group_by(Variable) %>%
  mutate(tmp = mean(tmp)) %>%
  arrange(Characteristic,.by_group = TRUE) %>%
  ungroup() %>%
  arrange(tmp) %>%
  select(-tmp) %>%
  mutate(across(Characteristic, 
                ~if_else(.x == "Immunosuppression",
                         "Immunosupp- ression",
                         .x))) %>%
  mutate(across(Variable, ~if_else(is.na(.x), " ", .x)))


table1_supplement <- list(
  col_names = col_names,
  cols_subtype = cols_subtype,
  data = table1_data_supplement
)

readr::write_rds(
  table1_supplement,
  here::here(release_folder, "table1_supplement.rds"))
