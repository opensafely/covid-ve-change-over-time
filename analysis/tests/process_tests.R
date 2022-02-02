################################################################################
# process tests data
library(tidyverse)
library(lubridate)
library(glue)

################################################################################
## source functions
# source(here::here("analysis", "lib", "data_properties.R"))

################################################################################
## create folders for outputs
fs::dir_create(here::here("output", "tests", "images"))
fs::dir_create(here::here("output", "tests", "tables"))

################################################################################
cat("--- read input_tests.feather ----")
data_tests_0 <- arrow::read_feather(
  file = here::here("output", "input_tests.feather")) %>%
  select(-starts_with(c("preg_36", "pregdel"))) %>%
  mutate(across(patient_id, as.integer))

data_eligible_e <- readr::read_csv(
  here::here("output", "data", "data_eligible_e.csv"))

################################################################################
# process pregnancy data
cat("--- process input_tests.feather ----")
data_pregnancy <- data_tests_0 %>%
  select(patient_id, starts_with("preg")) %>%
  pivot_longer(
    cols = -patient_id,
    names_patter = "preg_(\\d)",
    names_to = "comparison",
    values_to = "pregnancy"
  ) %>%
  mutate(across(c(comparison), factor))  

readr::write_rds(
  data_pregnancy,
  here::here("output", "data", "data_pregnancy.rds"),
  compress = "gz"
)

cat("--- process data_tests_1 ----")
data_tests_1 <- data_tests_0 %>%
  select(-starts_with("preg")) %>%
  mutate(across(contains("_date"), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days"))) 

################################################################################
# tabulate all vars 
# this is causing R to abort session- investigate
# data_properties(
#   data = data_tests,
#   path = file.path("output", "tests", "tables")
# )

################################################################################
# plot distibution of coviariates
cat("--- plot covariates ----")
plot_data <- data_tests_1 %>%
  transmute(patient_id,
         test_hist_1_n, 
         test_hist_2_n = test_hist_1_n + test_hist_2_n,
         test_hist_3_n) %>%
  pivot_longer(cols = -patient_id) %>%
  left_join(data_eligible_e, by = "patient_id") %>%
  mutate(across(arm, 
                factor, 
                levels = c("unvax", "vax"), 
                labels = c("unvaccinated", "vaccinated"))) %>%
  mutate(across(name,
                factor,
                labels = str_wrap(c("pre min(elig_date) in subtype", "pre elig_date", "elig_date to elig_date + 6 wks"), 25)))

# min for y axis 
min_y <- plot_data %>% 
  group_by(arm) %>%
  count() %>%
  ungroup() %>%
  mutate(p = 5/n) %>%
  summarise(p = max(p))

x_trunc <- 10

ggplot(NULL, aes(x = value)) +
  geom_bar(data = plot_data %>% filter(arm == "vaccinated"),
           aes(fill = arm, y = ..count../sum(..count..)), alpha = 0.5, width = 1) +
    geom_bar(data = plot_data %>% filter(arm == "unvaccinated"),
             aes(fill = arm, y = ..count../sum(..count..)), alpha = 0.5, width = 1) +
  scale_x_continuous(breaks = seq(1,10,1)) +
  scale_y_continuous(labels=scales::percent) +
  facet_wrap(~ name) +
  labs(y = "percent", x = "number of SARS-CoV-2 tests",
       caption = str_c(glue("x-axis truncated at {x_trunc}, y-axis truncated at "), signif(100*min_y$p,1), "% to mask bars corresponding to < 5 individuals")) +
  coord_cartesian(xlim = c(0,x_trunc), ylim = c(min_y$p, NA)) +
  scale_fill_discrete(name=NULL) +
  theme(legend.position = "bottom")
  
cat("--- save plot ----")
ggsave(
  filename = here::here("output", "tests", "images", "covariate_distribution.png"),
  width=15, height=20, units="cm")

################################################################################
# bin based on dist

data_tests_2 <- data_tests_1 %>%
  mutate(across(starts_with("test_hist"),
                ~ factor(case_when(
                  is.na(.x) ~ NA_character_,
                  .x < 1 ~ "0",
                  .x < 2 ~ "1",
                  .x < 3 ~ "2",
                  TRUE ~ "3+"
                )))) %>%
  select(-elig_date)

cat("--- save data_tests.rds ----")
readr::write_rds(
  data_tests_2,
  here::here("output", "data", "data_tests.rds"),
  compress = "gz"
)
