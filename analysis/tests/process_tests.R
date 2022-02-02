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

data_eligible_e <- readr::read_csv(
  here::here("output", "data", "data_eligible_e.csv"))

cat("--- process input_tests.feather ----")
data_tests <- data_tests_0 %>%
  mutate(across(contains("_date"), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days"))) #%>%
  # mutate(
  #   covid_test_both_elig_n = covid_test_pre_elig_n + covid_test_post_elig_n
  #   ) %>%
  # mutate(across(starts_with("covid_test"),
  #               ~ cut(.x, 
  #                     breaks = c(-Inf, 0, 4, Inf),
  #                     labels = c("0", "1-4", "5+"),
  #                     right = TRUE,
  #                     include.lowest = TRUE))) %>%
  # select(-elig_date)

# cat("--- check categorised variables ----")
# data_tests %>% select(starts_with("covid_test")) %>% summary()

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
cat("--- plot covariates ----")
plot_data <- data_tests %>%
  select(patient_id,
         covid_test_pre_elig_n,
         covid_test_post_elig_n) %>%
  pivot_longer(cols = -patient_id) %>%
  left_join(data_eligible_e, by = "patient_id") %>%
  mutate(across(arm, 
                factor, 
                levels = c("unvax", "vax"), 
                labels = c("unvaccinated", "vaccinated"))) %>%
  mutate(across(name,
                factor,
                levels = c("covid_test_pre_elig_n",
                           "covid_test_post_elig_n"),
                labels = c("pre 1st dose eligibility",
                           "6 weeks post 1st dose eligibility")))

# min for y axis 
min_y <- plot_data %>% 
  group_by(arm) %>%
  count() %>%
  ungroup() %>%
  mutate(p = 5/n) %>%
  summarise(p = max(p))

ggplot(NULL, aes(x = value)) +
  geom_bar(data = plot_data %>% filter(arm == "vaccinated"),
           aes(fill = arm, y = ..count../sum(..count..)), alpha = 0.5, width = 1) +
    geom_bar(data = plot_data %>% filter(arm == "unvaccinated"),
             aes(fill = arm, y = ..count../sum(..count..)), alpha = 0.5, width = 1) +
  scale_y_continuous(labels=scales::percent) +
  facet_wrap(~ name, scales = "free") +
  labs(y = "percent", x = "number of SARS-CoV-2 tests",
       caption = str_c("x-axis truncated at 20, y-axis truncated at ", signif(100*min_y$p,3), "% to mask bars corresponding to < 5 individuals")) +
  coord_cartesian(xlim = c(0,20), ylim = c(min_y$p, NA)) +
  scale_fill_discrete(name=NULL) +
  theme(legend.position = "bottom")
  
cat("--- save plot ----")
ggsave(
  filename = here::here("output", "tests", "images", "covariate_distribution.png"),
  width=15, height=20, units="cm")
