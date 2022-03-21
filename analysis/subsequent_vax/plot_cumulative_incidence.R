library(tidyverse)
library(glue)
library(survival)
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
  select(patient_id, covid_vax_2_date, start_of_period, arm)

# individuals eligible based on box d criteria
data_eligible_e_unvax <- readr::read_rds(
  here::here("output", "data", "data_eligible_e_unvax.rds"))  %>%
  mutate(arm = "unvax") %>%
  select(patient_id, start_of_period, arm)

# processed data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) %>%
  select(patient_id, subgroup, death_date, dereg_date)

data_wide_vax_dates <- readRDS(
  here::here("output", "data", "data_wide_vax_dates.rds")) %>%
  select(patient_id, covid_vax_1_date, covid_vax_3_date)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)
subgroup_order <- c(4,1,3,2)
subgroups_long <- if_else(
  subgroups %in% subgroups[c(2,3)],
  as.character(glue("{subgroups}*")),
  subgroups
)
subgroups_long_wrap <- str_wrap(subgroups_long, width = 25)

# redaction function for KM curves
source(here::here("analysis", "lib", "round_km.R"))

################################################################################

# if running locally read extracted data:
if(Sys.getenv("OPENSAFELY_BACKEND") %in% "") {
  
  release_folder <- "release20220226"
  image_path <- here::here(release_folder)
  
  survtable_redacted <- readr::read_csv(
    here::here(release_folder, "survtable_redacted.csv")) %>%
    mutate(across(subgroup,
                  ~ case_when(
                    str_detect(.x, "16-64") ~ 1,
                    str_detect(.x, "18-39") ~ 2,
                    str_detect(.x, "40-64") ~ 3,
                    str_detect(.x, "65") ~ 4,
                    TRUE ~ NA_real_
                  ))) %>%
    mutate(across(subgroup,
                  factor,
                  levels = subgroup_order,
                  labels = subgroups_long_wrap[subgroup_order])) 
  
} else { # else derive the data
  
  image_path <- here::here("output", "subsequent_vax", "images")
  
  # function to be applied in dplyr::filter
  no_evidence_of <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
  data_tte <- bind_rows(data_eligible_e_vax, data_eligible_e_unvax) %>%
    left_join(data_processed,
              by = "patient_id") %>%
    left_join(data_wide_vax_dates,
              by = "patient_id") %>%
    mutate(across(subgroup,
                  factor,
                  levels = subgroups[subgroup_order],
                  labels = subgroups_long_wrap[subgroup_order])) %>%
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
    # remove if subsequent vaccine, death or dereg before start_of_period
    filter_at(
      all_of(c("vax_date", "death_date", "dereg_date")),
      all_vars(no_evidence_of(., start_fu_date))) %>%
    group_by(start_fu_date) %>%
    mutate(across(c(end_fu_date, vax_date, death_date, dereg_date),
                  # time in weeks between start_fu_date and event
                  ~ as.integer(.x - start_fu_date)/7)) %>%
    ungroup() %>%
    select(-start_fu_date) %>%
    rename_at(vars(ends_with("_date")),
              ~ str_remove(.x, "_date")) %>%
    mutate(
      tte = pmin(vax, dereg, death, end_fu, na.rm = TRUE),
      status = if_else(
        !is.na(vax) & vax == tte,
        TRUE,
        FALSE
      )) %>%
    select(patient_id, arm, subgroup, tte, status) %>%
    mutate(across(arm, 
                  ~if_else(
                    arm == "unvax",
                    "Unvaccinated",
                    .x)
    ))
  
  survtable_redacted <- round_km(
    data = data_tte,
    time = "tte",
    event = "status",
    strata = c("subgroup", "arm")
  ) %>%
    ungroup() %>%
    mutate(
      c.inc = 1-surv,
      c.inc.ll = 1-surv.ul,
      c.inc.ul = 1-surv.ll
    ) %>%
    mutate(across(starts_with("c."), round, digits = 5)) %>%
    select(subgroup, arm, time, n.risk, n.event, n.censor, c.inc) 
  
  readr::write_csv(
    survtable_redacted,
    here::here("output", "subsequent_vax", "tables", "survtable_redacted.csv"))
  
}

################################################################################
# scale for x-axis
K <- study_parameters$max_comparisons
x_breaks <- seq(0,K*4,4)
x_labels <- x_breaks + 2

# create plot
plot_out <- survtable_redacted %>%
  mutate(across(arm, ~str_replace(.x, "ChAdOx", "ChAdOx1"))) %>%
  ggplot(aes(x = time, y = c.inc, colour = subgroup)) +
  geom_step(size=0.8) +
  facet_wrap(~ arm, nrow=2) +
  scale_color_viridis_d(
    name = "Subgroup"
  ) +
  scale_x_continuous(
    breaks = seq(0,24,4), # scale is time since start of period 1
    labels = seq(2,26,4) # label scale as time since second vax
  ) +
  labs(
    x = "Weeks since second dose",
    y = "Incidence of 1st dose                 Incidence of 3rd dose"
    ) +
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
    
    legend.position = c(0.75,0.2),
    legend.key.size = unit(0.8, "cm")
    
  )

ggsave(plot = plot_out,
       filename = file.path(image_path, "ci_vax.png"),
       width=16, height=12, units="cm")
  