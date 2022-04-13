# rmarkdown:::with_pandoc_safe_environment(
#   system(paste(shQuote(Sys.which("pandoc")), "--version"), intern = TRUE)
# )

# Sys.which("pandoc")
cat("------check pandoc")
rmarkdown::find_pandoc()

library(tidyverse)
library(lubridate)
library(kableExtra)

################################################################################
# create output folder
fs::dir_create(here::here("output", "eda"))

# rmarkdown::render(
#   "analysis/eda/covidadmitted_postest.Rmd",
#   output_file="covidadmitted_postest",
#   output_dir="output/eda")

################################################################################
# read data for ever covariates
data_ever <- arrow::read_feather(
  file = here::here("output", "input_ever.feather")) 

data_0 <- data_ever %>%
  select(patient_id, covidadmitted_post_date, starts_with("covidadmitted_postest")) %>%
  mutate(across(contains("_date"), 
                ~ floor_date(
                  as.Date(.x, format="%Y-%m-%d"),
                  unit = "days"))) 

################################################################################
capture.output(
  {
    data_0 %>%
      mutate(covidadmitted = !is.na(covidadmitted_post_date)) %>%
      group_by(covidadmitted) %>%
      count() %>%
      ungroup() %>%
      mutate(percent = 100*n/sum(n)) %>%
      kable(
        "pipe",
        caption = "Percent of patients with a COVID-19 hospital admission")
  },
  file = here::here("output", "eda", "covidadmitted_postest.txt"),
  append = FALSE
)

################################################################################
data_1 <- data_0 %>%
  filter(!is.na(covidadmitted_post_date)) %>%
  pivot_longer(
    cols = starts_with("covidadmitted_postest"),
    names_pattern = "covidadmitted_postest_(.)_date",
    names_transform = as.integer,
    values_drop_na = TRUE
  ) %>%
  mutate(diff = as.numeric(value - covidadmitted_post_date)) %>%
  group_by(patient_id) %>%
  mutate(
    min_postest_date = min(name),
    max_postest_date = max(name)
    ) %>%
  ungroup()

################################################################################  
# percent of patients with i=1-6 positive tests
capture.output(
  {
    data_1 %>% 
      group_by(name) %>%
      count() %>%
      ungroup() %>%
      mutate(percent = 100*n/sum(n)) %>%
      kable(
        "pipe",
        caption = "Number of positive tests per patient during the window"
      )
  },
  file = here::here("output", "eda", "covidadmitted_postest.txt"),
  append = TRUE
)

# distribution of all positive test dates relative to date of hospital admission
data_1 %>%
  ggplot(aes(x = diff)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_histogram(
    center = 0.5,
    binwidth = 1, 
    alpha = 0.8) +
  scale_x_continuous(
    breaks = seq(-42,7,7),
    limits = c(-42,7),
    labels = as.character(seq(-6,1,1))
  ) +
  labs(
    x = "Weeks since COVID-19 hospital admission",
    y = "Number of positive tests"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(margin = margin(t=10)),
    axis.title.y = element_text(margin = margin(r=10))
  )
ggsave(
  filename = here::here("output", "eda", "covidadmitted_postest_all.png")
)

# distribution of min positive test dates relative to date of hospital admission
data_1 %>%
  filter(name == min_postest_date) %>%
  ggplot(aes(x = diff)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_histogram(
    center = 0.5,
    binwidth = 1, 
    alpha = 0.8) +
  scale_x_continuous(
    breaks = seq(-42,7,7),
    limits = c(-42,7),
    labels = as.character(seq(-6,1,1))
  ) +
  labs(
    x = "Weeks since COVID-19 hospital admission",
    y = "Number of positive tests"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(margin = margin(t=10)),
    axis.title.y = element_text(margin = margin(r=10))
  )
ggsave(
  filename = here::here("output", "eda", "covidadmitted_postest_min.png")
)