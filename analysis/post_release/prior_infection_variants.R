# check possible variants for each of the time since prior infection categories

################################################################################
library(tidyverse)
library(glue)
library(flextable)
library(officer)

################################################################################
# release folder
if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")  

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

################################################################################
# read data
check_fu <- bind_rows(
  lapply(
    subgroups,
    function(x)
      readr::read_delim(here::here(release_folder, "variants", glue("check_fu_{x}.txt")), 
                        delim = "|", escape_double = FALSE, trim_ws = TRUE) %>%
      select(k,date,n) %>%
      slice(-1) %>%
      mutate(subgroup = x) 
  )
) 

################################################################################
# clean data
data <- check_fu %>%
  mutate(across(subgroup, factor, levels = subgroups)) %>%
  mutate(across(c(k,n), as.integer)) %>%
  mutate(across(date, as.Date)) %>%
  filter(n>0) %>%
  group_by(subgroup, k) %>%
  summarise(start_date = min(date), end_date = max(date), .groups = "keep") %>%
  ungroup() %>%
  mutate(
    wild_start_date = as.Date("2020-03-01"),
    alpha_start_date = as.Date("2020-12-01"),
    delta_start_date = as.Date("2021-06-01"),
    # min and max time since variant
    # wild
    wild_min = as.integer(start_date - alpha_start_date - 1),
    wild_max = as.integer(end_date - wild_start_date),
    # alpha
    alpha_min = as.integer(start_date - delta_start_date - 1),
    alpha_max = as.integer(end_date - alpha_start_date),
    # delta
    delta_min = as.integer(start_date - as.Date(study_parameters$end_date)),
    delta_max = as.integer(end_date - delta_start_date)
  ) %>%
  select(-ends_with("date")) %>%
  mutate(across(starts_with(c("wild", "alpha", "delta")),
                ~case_when(
                  .x < 90 ~ "A",
                  .x < 180 ~ "B",
                  TRUE ~ "C"))) %>%
  pivot_longer(
    cols = starts_with(c("wild", "alpha", "delta")),
    names_to = c(".value", "type"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  select(-type) %>%
  distinct() %>%
  pivot_longer(
    cols = starts_with(c("wild", "alpha", "delta"))
  ) %>%
  distinct() %>%
  mutate(tmp = TRUE) %>%
  pivot_wider(
    names_from = c(name,value),
    values_from = tmp
  ) %>%
  mutate(across(where(is.logical), ~if_else(is.na(.x), FALSE, .x)))

# all possible categories
tmp <- bind_cols(
  unlist(
    lapply(
      c("wild", "alpha", "delta"), 
      function(x) 
        lapply(
          c("A", "B", "C"), 
          function(y) 
            tibble(!! as.character(glue("{x}_{y}")) := FALSE)
          )
      ), 
    recursive = FALSE)
  ) 

# prepare data for table
table <- data %>%
  bind_cols(tmp %>% select(-any_of(names(data)))) %>%
  select(subgroup, k, names(tmp)) %>%
  mutate(across(wild_B, ~if_else(wild_A & wild_C, TRUE, .x))) %>%
  mutate(across(alpha_B, ~if_else(alpha_A & alpha_C, TRUE, .x))) %>%
  mutate(across(delta_B, ~if_else(delta_A & delta_C, TRUE, .x))) %>%
  mutate(across(starts_with(c("wild", "alpha", "delta")), 
                ~if_else(.x, "yes", NA_character_)))

# key for table names
col_names <- names(table)
col_names_clean <- str_remove(col_names, "wild_|alpha_|delta_")
col_names_clean <- str_replace(col_names_clean, "A", "1-90")
col_names_clean <- str_replace(col_names_clean, "B", "91-180")
col_names_clean <- str_replace(col_names_clean, "C", "181+")
names(col_names_clean) <- col_names

# generate table
flextable1 <- table %>%
  flextable() %>%
  set_header_labels(
    values = as.list(col_names_clean)
  ) %>%
  merge_v(j=1, part = "body") %>%
  add_header_row(
    values = c("", "wild", "alpha", "delta"),
    colwidths = c(2, rep(3,3))
  ) %>%
  hline(
    i = seq(6,4*6,6)
  ) %>%
  vline(
    j = seq(2, 3*3+2,3),
    part = "all"
  ) %>%
  border_outer(part = "all")

# save table as docx
doc <- read_docx() %>%
  body_add_flextable(value = flextable1, split = FALSE) %>%
  print(target = here::here(release_folder, "prior_infection_variants_table.docx"))

