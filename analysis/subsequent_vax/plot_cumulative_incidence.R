library(tidyverse)
library(RColorBrewer)
library(glue)
library(survival)
library(survminer)
library(lubridate)

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
  select(patient_id, start_of_period, arm)

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
    start_fu_date = start_of_period + days(14),
   # end date of final comparison or end of data availability
   end_fu_date = pmin(start_fu_date + days(study_parameters$max_comparisons*28),
                      study_parameters$end_date)
    ) %>%
  select(-start_of_period) %>%
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
    select(patient_id, arm, subgroup, tte, status)
  
}

################################################################################


plot_fun <- function(data_tte, censor_var) {
  
  fit <- survfit(Surv(tte, status) ~ arm + subgroup, 
                 data = data_tte)
  
  # Plot cumulative events
  survplots <- ggsurvplot(fit,
                          data = data_tte,
                          break.time.by = 4,
                          xlim = c(0, max(data_tte$tte)),
                          conf.int = FALSE,
                          censor=FALSE, # show censor ticks on line?
                          cumevents = TRUE,
                          cumcensor = TRUE,
                          fun = "event",
                          linetype = "arm",
                          size = 0.5,
                          color = "subgroup",
                          palette = brewer.pal(length(subgroups), name = "Set2"),
                          legend.title = "",
                          font.x = 10,
                          font.y = 10,
                          font.tickslab = 10,
                          font.title = 12,
                          font.caption = 10,
                          font.legenf = 10)

  if (censor_var == "postest") {
    censor_var_long <- "positive COVID-19 test, "
  } else if (censor_var == "covidadmitted") {
    censor_var_long <- "COVID-19 hospital admission, "
  } else {
    censor_var_long <- ""
  }

  caption_string <- glue("Events are 1st dose in unvaccinated arm and 3rd dose in BNT162b2 and ChAdOx arms. Follow-up is censored at {censor_var_long}death and de-registration.")

  plot_out <- survplots$plot +
    labs(
      title = "Cumulative incidence of 1st and 3rd vaccine doses",
      x = "weeks since start of comparison period 1",
      y = "cumulative event",
      caption = str_wrap(caption_string, 80)
    ) +
    scale_linetype_discrete(name = NULL) +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          plot.caption = element_text(hjust = 0))

  ggsave(plot = plot_out,
         filename = here::here("output", "subsequent_vax", "images", glue("ci_vax_{censor_var}.png")),
         width=15, height=12, units="cm")

  survtable <- survplots$data.survtable %>%
    select(arm, subgroup, time, n.risk, cum.n.event, cum.n.censor) %>%
    # round to the nearest 10
    mutate(across(c(n.risk, cum.n.event, cum.n.censor),
                  ~round(.x,-1)))

  capture.output(
    survtable %>%
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

