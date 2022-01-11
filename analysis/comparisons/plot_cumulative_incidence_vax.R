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

data_vax_incidence <- bind_rows(
  readr::read_rds(
    here::here("output", "data", glue("data_comparisons_{brand}_vax.rds"))) %>%
    mutate(arm = "vax"),
  readr::read_rds(
    here::here("output", "data", glue("data_comparisons_{brand}_unvax.rds"))) %>%
    mutate(arm = "unvax")
) %>%
  mutate(subgroup = case_when(
    jcvi_group %in% "02" ~ "02",
    jcvi_group %in% c("11", "12") ~ "11-12",
    TRUE ~ "03-10"
  )) %>%
  filter(subgroup %in% subgroups) %>%
  group_by(patient_id, subgroup, arm) %>%
  summarise(
    min_start_fu_date = min(start_fu_date),
    max_end_fu_date = max(end_fu_date),
    .groups = "keep"
    ) %>%
  ungroup() %>%
  left_join(
    readr::read_rds(
      here::here("output", "data", "data_wide_vax_dates.rds")
    ),
    by = "patient_id"
  ) %>%
  left_join(
    readr::read_rds(
      here::here("output", "data", "data_covid_any.rds")),
    by = "patient_id"
  ) %>%
  mutate(
    vax_date = if_else(
      arm %in% "vax",
      covid_vax_3_date,
      covid_vax_1_date
    ),
    end_at = pmin(
      vax_date, # date of 1st of 3rd dose
      covid_any_date, # date of any covid event
      max_end_fu_date, # end of individual's followup
      as.Date(study_parameters$end_date), # last date of data
      na.rm = TRUE),
    time = as.integer(end_at - min_start_fu_date),
    status = case_when(
      is.na(vax_date) ~ 0L,
      vax_date > as.Date(end_at) ~ 0L,
      TRUE ~ 1L
    )
  ) %>%
  mutate(across(arm, ~if_else(.x %in% "vax", "vaccinated", "unvaccinated"))) %>%
  mutate(strata = str_c(subgroup, "; ", arm)) %>%
  select(patient_id, strata, time, status)

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
  
