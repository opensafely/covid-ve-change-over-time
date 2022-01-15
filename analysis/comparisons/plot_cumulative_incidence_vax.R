library(tidyverse)
library(RColorBrewer)
library(glue)
library(survival)
library(survminer)

## import study_parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  brand <- "BNT162b2"
  
} else{
  brand <- args[[1]]
}

fs::dir_create(here::here("output", "images"))

################################################################################

if (brand %in% "BNT162b2") {
  subgroups <- c("02", "03-10", "11-12")
} else {
  subgroups <- "03-10"
}

################################################################################

data_comparisons <- bind_rows(
  readr::read_rds(
    here::here("output", "comparisons", "data", "data_comparisons_BNT162b2.rds")),
  readr::read_rds(
    here::here("output", "comparisons", "data", "data_comparisons_ChAdOx.rds")),
  readr::read_rds(
    here::here("output", "comparisons", "data", "data_comparisons_unvax.rds"))
) 

data_wide_vax_dates <- readRDS(
  here::here("output", "data", "data_wide_vax_dates.rds")) %>%
  select(patient_id, covid_vax_1_date, covid_vax_3_date)

data_min <- data_comparisons %>%
  # only keep earliest date for each individual
  group_by(patient_id) %>%
  summarise(across(
    c(start_fu_date,
      postest_date,covidadmitted_date,
      death_date, dereg_date),
    min, na.rm = TRUE),
    .groups = "keep") %>%
  ungroup() 

data_max <- data_comparisons %>%
  # only keep earliest date for each individual
  group_by(patient_id) %>%
  summarise(across(
    end_fu_date,
    max, na.rm = TRUE),
    .groups = "keep") %>%
  ungroup() 

data_tte <- data_comparisons %>%
  distinct(patient_id, arm, subgroup) %>%
  left_join(data_min, by = "patient_id") %>%
  left_join(data_max, by = "patient_id") %>%
  left_join(data_wide_vax_dates, by = "patient_id") %>%
  mutate(vax_date = if_else(
    arm %in% "unvax",
    covid_vax_1_date,
    covid_vax_3_date)) %>%
  select(-covid_vax_1_date, -covid_vax_3_date) %>%
  group_by(arm, subgroup) %>%
  mutate(across(c(start_fu_date, end_fu_date, vax_date,
                  postest_date,covidadmitted_date,
                  death_date, dereg_date),
                ~ as.integer(.x - min(start_fu_date)))) %>%
  ungroup() %>%
  rename_at(vars(ends_with("_date")),
            ~ str_remove(.x, "_date"))


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
    select(patient_id, arm, subgroup, tstart = start_fu, tstop = tte, status)
  
}

data_tte %>% censor_fun()
data_tte %>% censor_fun(postest)


################################################################################

fit <- survfit(Surv(time, status) ~ strata, 
               data = data_vax_incidence)

################################################################################
  
caption_width <- 100
exdent_width <- 10

# Plot cumulative events
survplots <- ggsurvplot(fit, 
                        break.time.by = 4,
                        xlim = c(0, max(data_vax_incidence$time)),
                        conf.int = TRUE,
                        # palette = brewer.pal(n = length(subgroups), name = "Dark2"),
                        censor=TRUE, # show censor ticks on line?
                        cumevents = TRUE, 
                        cumcensor = TRUE, 
                        # risk.table.col = "strata",
                        fun = "event",
                        # aesthetics
                        break.x.by = 28,
                        xlab = "days since start of individual's follow-up",
                        legend.title = "JCVI group(s); arm",
                        legend = "bottom",
                        font.legend = 8,
                        # legend.labs = strata,
                        caption = str_c(
                          str_wrap(
                            "For the unvaccinated arm: start of follow-up is 14 days after the start of the second vaccination period for the individual's strata, and an event is a first dose of vaccination.",
                            caption_width, 
                            exdent = exdent_width),
                          "\n",
                          str_wrap(
                            "For the vaccinated arm: start of follow-up is 14 days after the individual's date of second vaccination, and an event is a third dose of vaccination.",
                            caption_width, 
                            exdent = exdent_width)
                          ),
                        ggtheme = theme_bw() + 
                          theme(plot.caption = element_text(hjust = 0, size=8)))

################################################################################

ggsave(plot = survplots$plot,
       filename = here::here("output",  "images", glue("cumulative_incidence_{brand}.png")),
       width=15, height=12, units="cm")

################################################################################

survtable <- survplots$data.survtable %>%
  select(strata, time, n.risk, cum.n.event, cum.n.censor) %>%
  # round to the nearest 10
  mutate(across(c(n.risk, cum.n.event, cum.n.censor), 
                ~round(.x,-1))) 

capture.output(
  survtable %>%
    kableExtra::kable("pipe", padding = 2),
  file = here::here("output", "tables", glue("survtable_{brand}.txt")),
  append=FALSE
)

################################################################################
  
