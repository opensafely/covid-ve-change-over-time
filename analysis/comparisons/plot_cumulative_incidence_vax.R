

library(tidyverse)
library(RColorBrewer)
library(glue)
library(survival)
library(survminer)

## import study_parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

comparison <- "BNT162b2"
arm <- "vax"

if (arm == "vax") censor_at <- "covid_vax_3_date" else censor_at <- "covid_vax_1_date" 


data_vax_incidence <- readr::read_rds(
  here::here("output", "data", glue("data_comparisons_{comparison}_{arm}.rds"))) %>%
  mutate(subgroup = case_when(
    jcvi_group %in% "02" ~ "02",
    jcvi_group %in% c("11", "12") ~ "11-12",
    TRUE ~ "03-10"
  )) %>%
  select(patient_id, arm, subgroup, start_fu_date) %>%
  # keep earliest start_fu_date for each individual
  arrange(patient_id, start_fu_date) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  left_join(
    readr::read_rds(
      here::here("output", "data", "data_wide_vax_dates.rds")
    ),
    by = "patient_id"
  ) %>%
  mutate(
    censor_at = pmin(!! sym(censor_at), as.Date(study_parameters$end_date)),
    time = as.integer(censor_at - start_fu_date),
    status = as.integer(!is.na(!! sym(censor_at)))
    ) %>%
  select(patient_id, arm, subgroup, time, status)



cat("#### fit survival model ####\n")
fit <- survfit(Surv(time, status) ~ subgroup, 
               data = data_vax_incidence)

# write_rds(fit, here::here("output", "models", glue("surv_model_{group}.rds")))

subgroups <- unique(data_vax_incidence$subgroup)

cat("#### generate plots ####\n")
# Plot cumulative events
survplots <- ggsurvplot(fit, 
                        break.time.by = 4,
                        xlim = c(0,max(data_vax_incidence$time)),
                        conf.int = TRUE,
                        palette = brewer.pal(n = length(subgroups), name = "Dark2"),
                        censor=TRUE, #don't show censor ticks on line
                        cumevents = FALSE, 
                        cumcensor = FALSE, 
                        risk.table.col = "strata",
                        fun = "event",
                        # aesthetics
                        break.x.by = 28,
                        xlab = "Time since second dose + 14 days",
                        # legend.title = title,
                        # legend.labs = strata,
                        ggtheme = theme_bw())
  



  
  
  
  group_by(subgroup, jcvi_group, elig_date, region) %>%
  summarise(start_fu_date = min(start_fu_date))
