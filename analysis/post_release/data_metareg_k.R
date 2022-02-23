library(tidyverse)

metareg_results_0 <- readxl::read_excel(here::here("release20220221", "metareg_results.xlsx"))

metareg_results <- metareg_results_0 %>%
  mutate(across(outcome,
                ~case_when(.x %in% "Any test" ~ "anytest",
                           .x %in% "Positive test" ~ "postest",
                           .x %in% "COVID-19 hospitalisation" ~ "covidadmitted",
                           .x %in% "COVID-19 death" ~ "coviddeath",
                           .x %in% "Non-COVID death" ~ "noncoviddeath",
                           TRUE ~ NA_character_))) %>%
  mutate(
    # slope ci
    logrhr_lower = logrhr - qnorm(0.975)*selogrhr,
    logrhr_higher = logrhr + qnorm(0.975)*selogrhr,
    # intercept ci
    loghr1_lower = loghr1 - qnorm(0.975)*seloghr1,
    loghr1_higher = loghr1 + qnorm(0.975)*seloghr1,
    ) 

metareg_results_k <- metareg_results %>%
  distinct(subgroup, comparison, outcome) %>%
  mutate(k=factor(1, levels=1:6)) %>%
  complete(subgroup, comparison, outcome, k) %>%
  left_join(
    metareg_results,
    by = c("subgroup", "comparison", "outcome")
  ) %>%
  mutate(across(k, ~as.integer(as.character(.x)))) %>%
  mutate(
    line = loghr1 + (k-1)*logrhr#, # use k-1 because intercept at k=-1
    # CI for slope not valid as don't know covariance between intercept and slope
    # line_lower = loghr1_lower + (k-1)*logrhr_lower,
    # line_higher = loghr1_higher + (k-1)*logrhr_higher
  ) %>%
  mutate(across(starts_with("line"), exp))

readr::write_rds(
  metareg_results_k,
  here::here("release20220221", "metareg_results_k.rds")
)

# metareg_results_k %>%
#   filter(subgroup =="65+ years", comparison!="both", outcome == "COVID-19 hospitalisation") %>%
#   ggplot(aes(x = k)) +
#   geom_line(aes(y = line, colour=comparison)) +
#   geom_ribbon(aes(ymin=line_lower, ymax=line_higher, fill = comparison), alpha=0.2)
