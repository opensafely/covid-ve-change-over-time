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
  
  data_tte_list <- data_tte %>% arrange(arm) %>% group_split(arm)
  
  names_data_tte_list <- sapply(
    data_tte_list,
    function(x) unique(x$arm)
  )
  names_data_tte_list[which(names_data_tte_list=="unvax")] <- "unvaccinated"
  
  apply_ggsurvplot <- function(.data) {
    
    fit <- survfit(Surv(tte, status) ~ arm + subgroup, 
                   data = .data)
    
    ggsurvplot(fit,
               data = .data,
               break.time.by = 4,
               # xlim = c(0, max(.data$tte)),
               conf.int = FALSE,
               censor=FALSE, # show censor ticks on line?
               cumevents = TRUE,
               cumcensor = TRUE,
               fun = "event",
               size = 0.5,
               color = "subgroup",
               palette = brewer.pal(length(subgroups), name = "Set2"),
               legend.title = "Subgroup",
               font.x = 8,
               font.y = 8,
               font.tickslab = 10,
               font.title = 12,
               font.caption = 10,
               font.legenf = 10)
  }
  
  # scale for x-axis
  K <- study_parameters$max_comparisons
  x_breaks <- seq(0,K*4,4)
  x_labels <- x_breaks + 2
  
  survplots <- lapply(data_tte_list, apply_ggsurvplot)
  
  max_ci_strata <- bind_rows(
    lapply(
      survplots,
      function(x)
        x$data.survplot %>% 
        filter(time == max(time)) %>%
        mutate(ci = 1-surv) %>%
        select(strata, ci)
      )
    )
  max_ci <- max(max_ci_strata$ci)
  
  plot_legend <- get_legend(
    survplots[[1]]$plot +
      theme(
        legend.direction = "vertical"
        )
    )
  
  survplots_processed <- lapply(
    survplots,
    function(x)
      x$plot + 
      scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels
        ) +
      lims(y = c(0,max_ci)) +
      theme(
        plot.title = element_text(size=10),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        legend.position = "none"
        )
  )
  
  plot_out <- plot_grid(survplots_processed[[1]] + 
              labs(title = glue("3rd dose in {names_data_tte_list[1]} arm"),
                   x = "weeks since second dose\n"), 
            survplots_processed[[2]] +
              labs(title = glue("3rd dose in {names_data_tte_list[2]} arm"),
                   x = "weeks since second dose\n"), 
            survplots_processed[[3]] + 
              labs(title = glue("1st dose in {names_data_tte_list[1]} arm"),
                   x = "weeks since start of\nsecond vaccination period"),
            plot_legend,
            labels = c("A", "B", "C"))
  
  if (censor_var == "postest") {
    censor_var_long <- "positive COVID-19 test, "
  } else if (censor_var == "covidadmitted") {
    censor_var_long <- "COVID-19 hospital admission, "
  } else {
    censor_var_long <- ""
  }

  caption_string <- glue("Follow-up is censored at {censor_var_long}death and de-registration.")
  
  plot_out <- ggdraw(add_sub(plot_out, caption_string, x = 0, hjust = 0, size=10))

  ggsave(plot = plot_out,
         filename = here::here("output", "subsequent_vax", "images", glue("ci_vax_{censor_var}.png")),
         width=15, height=12, units="cm")
  
  survtable <- bind_rows(
    lapply(
      survplots,
      function(x)
        x$data.survtable %>% 
        mutate(less5 = if_else(
          n.event > 0 & n.event <=5,
          TRUE,
          FALSE
        )) %>%
        select(arm, subgroup, time, n.risk, n.event, cum.n.event, cum.n.censor, less5)
        
    )
  )


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

