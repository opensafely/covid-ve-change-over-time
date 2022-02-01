################################################################################
# process tests data
library(tidyverse)
library(lubridate)

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
  file = here::here("output", "input_tests.feather")) 

cat("--- process input_tests.feather ----")
data_tests <- data_tests_0 %>%
  mutate(across(contains("_date"), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days"))) %>%
  mutate(
    covid_test_both_elig_n = covid_test_pre_elig_n + covid_test_post_elig_n
    ) %>%
  mutate(across(starts_with("covid_test"),
                ~ cut(.x, 
                      breaks = c(-Inf, 0, 4, Inf),
                      labels = c("0", "1-4", "5+"),
                      right = TRUE,
                      include.lowest = TRUE))) %>%
  select(-elig_date)

cat("--- check categorised variables ----")
data_tests %>% select(starts_with("covid_test")) %>% summary()


cat("--- save data_tests.rds ----")
readr::write_rds(
  data_tests,
  here::here("output", "data", "data_tests.rds"),
  compress = "gz"
)

################################################################################
# tabulate all vars 
# this is causing R to abort session- investigate
# data_properties(
#   data = data_tests,
#   path = file.path("output", "tests", "tables")
# )

################################################################################
# plot distibution of coviariates
# cat("--- plot covariates ----")
# data_tests %>%
#   select(patient_id, 
#          pre=covid_test_pre_elig_n,
#          post=covid_test_post_elig_n) %>%
#   mutate(both = pre + post) %>%
#   pivot_longer(cols = -patient_id) %>%
#   mutate(across(value, ~if_else(.x > 30, 30L, .x))) %>%
#   ggplot(aes(x = value)) +
#   geom_bar() +
#   facet_wrap(~ name, scales = "free") +
#   scale_y_log10() +
#   theme(legend.position = "bottom")
# cat("--- save plot ----")
# ggsave(
#   filename = here::here("output", "tests", "images", "covariate_distribution.png"),
#   width=20, height=15, units="cm")