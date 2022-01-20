library(tidyverse)
library(RColorBrewer)
library(glue)
library(survival)
library(survminer)
library(lubridate)
library(cowplot)

################################################################################

fs::dir_create(here::here("output", "subsequent_vax", "images"))
fs::dir_create(here::here("output", "subsequent_vax", "tables"))

## import study_parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))


# individuals eligible based on box c criteria
data_eligible_e_vax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_vax.rds"))  %>%
  mutate(arm = brand) %>%
  select(patient_id, covid_vax_2_date, start_of_period, arm)

# individuals eligible based on box d criteria
data_eligible_e_unvax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_unvax.rds"))  %>%
  mutate(arm = "unvax") %>%
  select(patient_id, start_of_period, arm)

# processed data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) %>%
  select(patient_id, subgroup,
         postest_date, covidadmitted_date, death_date, dereg_date)

data_wide_vax_dates <- readRDS(
  here::here("output", "data", "data_wide_vax_dates.rds")) %>%
  select(patient_id, covid_vax_1_date, covid_vax_3_date)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))

# redaction functions
source(here::here("analysis", "lib", "redaction_functions.R"))

source(here::here("analysis", "lib", "data_process_functions.R"))

################################################################################

# function to be applied in dplyr::filter
no_evidence_of <- function(cov_date, index_date) {
  is.na(cov_date) | index_date < cov_date
}

censor_vars <- c("death_date", "dereg_date")

data <- bind_rows(data_eligible_e_vax, data_eligible_e_unvax) %>%
  left_join(data_processed,
            by = "patient_id") %>%
  left_join(data_wide_vax_dates,
            by = "patient_id") %>%
  mutate(
    # start date of comparison 1 
    start_fu_date = if_else(
      arm %in% "unvax",
      start_of_period + days(14),
      covid_vax_2_date + days(14)
      ),
   # end date of final comparison or end of data availability
   end_fu_date = pmin(start_fu_date + days(study_parameters$max_comparisons*28),
                      study_parameters$end_date)
    ) %>%
  select(-start_of_period, -covid_vax_2_date) %>%
  mutate(vax_date = if_else(
    arm %in% "unvax",
    covid_vax_1_date,
    covid_vax_3_date)) %>%
  select(-covid_vax_1_date, -covid_vax_3_date) %>%
  # remove if death or dereg before start_of_period
  filter_at(
    all_of(censor_vars),
    all_vars(no_evidence_of(., start_fu_date))) %>%
  group_by(start_fu_date) %>%
  mutate(across(c(end_fu_date, vax_date,
                  postest_date,covidadmitted_date,
                  death_date, dereg_date),
                # time in weeks between start_fu_date and event
                ~ as.integer(.x - start_fu_date)/7)) %>%
  ungroup() %>%
  select(-start_fu_date) %>%
  rename_at(vars(ends_with("_date")),
            ~ str_remove(.x, "_date"))

################################################################################

censor_fun <- function(.data, ...) {
  
  dots <- enquos(...) # get a list of quoted dots
  
  .data %>%
    mutate(
      tte = pmin(!!! dots,  # unquote-splice the list
                 vax, dereg, death, end_fu, na.rm = TRUE),
      status = if_else(
        !is.na(vax) & vax == tte,
        TRUE,
        FALSE
      )) %>%
    # postest and covid admitted may occur before zero
    # if censoring on postest, must remove individuals with postest before start of comparison 1
    # likewise for covidadmitted
    # when not censoring on either, tte will not be <=0, 
    # as neither postest nor covid admitted included in tte calculation in pmin above, 
    # so no individuals removed by tte>0
    filter(tte > 0) %>%
    select(patient_id, arm, subgroup, tte, status) %>%
    mutate(group = if_else(
      arm == "unvax",
      "unvaccinated",
      "vaccinated"
    ))
  
}

################################################################################


plot_fun <- function(data_tte, censor_var) {
  
  # scale for x-axis
  K <- study_parameters$max_comparisons
  x_breaks <- seq(0,K*4,4)
  x_labels <- x_breaks + 2
  
  fit <- survfit(Surv(tte, status) ~ arm + subgroup, 
                 data = data_tte)
  
  ggsurvplot_res <- ggsurvplot(fit,
                               data = data_tte,
                               break.time.by = 4)  
  
  if (censor_var == "postest") {
    censor_var_long <- "positive COVID-19 test, "
  } else if (censor_var == "covidadmitted") {
    censor_var_long <- "COVID-19 hospital admission, "
  } else {
    censor_var_long <- ""
  }
  
  caption_string <- glue("For the unvaccinated arm the time-scale is \"weeks since start of second vaccination period\". Follow-up for all arms is censored at {censor_var_long}death and de-registration.")
  
  
  survtable <- ggsurvplot_res$data.survplot %>%
    group_by(strata) %>%
    # suppress small numbers
    mutate(across(c(n.event, n.censor),
                  ~ case_when(
                    0 <= .x & .x < 3 ~ 0,
                    .x < 5 ~ 5,
                    TRUE ~ .x
                  ))) %>%
    mutate(n.excluded = lag(cumsum(n.event + n.censor)),
           n.risk_new = max(n.risk) - n.excluded) %>%
    # recalculate n.risk now that small numbers suppressed in n.event and n.censor
    mutate(across(n.risk,
                  ~ if_else(
                    is.na(n.risk_new),
                    .x,
                    n.risk_new)
    )) %>%
    ungroup() %>%
    select(arm, subgroup, time, n.risk, n.event, n.censor, surv)
  
  
  plot_out <- survtable %>%
    mutate(across(arm,
                  ~ fct_case_when(
                    .x %in% "BNT162b2" ~ "3rd dose in BNT162b2 arm",
                    .x %in% "ChAdOx" ~ "3rd dose in ChAdOx arm",
                    .x %in% "unvax" ~ "1st dose in unvaccinated arm"
                  )
                  )) %>%
    ggplot(aes(x = time, y = 1-surv, colour = subgroup)) +
    geom_step() +
    facet_wrap(~ arm, nrow = 2) +
    labs(x = "weeks since second dose",
         y = "cumulative event",
         caption = str_wrap(caption_string,102),
         title = "Cumulative incidence of subsequent vaccination") +
    scale_color_discrete(name = NULL) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=8),
      
      axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 10, l = 0)),
      axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0, size = 11),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = c(0.75,0.25)
      
    )
  
  ggsave(plot = plot_out,
         filename = here::here("output", "subsequent_vax", "images", glue("ci_vax_{censor_var}.png")),
         width=15, height=12, units="cm")
  

  capture.output(
    survtable %>%
      select(-surv) %>%
      kableExtra::kable("pipe", padding = 2),
    file = here::here("output", "subsequent_vax", "tables", glue("survtable_{censor_var}.txt")),
    append=FALSE
  )
  
}

################################################################################

data_tte_none <- data %>% censor_fun()
data_tte_postest <- data %>% censor_fun(postest)
data_tte_covidadmitted <- data %>% censor_fun(covidadmitted)

plot_fun(data_tte = data_tte_none, "none")
plot_fun(data_tte = data_tte_postest, "postest")
plot_fun(data_tte = data_tte_covidadmitted, "covidadmitted")

