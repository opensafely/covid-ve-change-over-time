################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobs <- FALSE
  group <- "02"
  outcome <- "postest"
  
} else{
  removeobs <- TRUE
  group <- args[[1]]
  outcome <- args[[2]]
}

################################################################################
# study_parameters <- readr::read_rds(
#   here::here("output", "lib", "study_parameters.rds"))


fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "tables"))
fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "images"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds")) %>%
  filter(jcvi_group %in% group, include) %>%
  distinct(brand, n_comparisons)


formatpercent100 <- function(x,accuracy){
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}


for (b in as.character(unique(second_vax_period_dates$brand))) {
  
  model_tidy <- readr::read_rds(
    here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{outcome}_modelcox_tidy.rds"))) 
  
  K <- second_vax_period_dates$n_comparisons[second_vax_period_dates$brand == b]
  
  ends <- seq(14, (K+1)*28, 28)
  starts <- ends + 1
  days_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")
  
  plot_data <- model_tidy %>% 
    filter(str_detect(term, "^comparison"), !str_detect(label, "unvax")) %>%
    mutate(
      comparison = factor(as.integer(str_extract(label, "\\d")),
                          labels = days_since_2nd_vax),
      
      # double check if estimates need exponentiated!
      # estimate = exp(estimate),
      # lower = exp(conf.low),
      # upper = exp(conf.high)
      
      estimate = estimate,
      lower = conf.low,
      upper = conf.high
      
    ) %>%
    mutate(across(model,
                  factor,
                  levels = c(0,1,2),
                  labels = c("Region-stratfied Cox model, no further adjustment", 
                             "Region-stratfied Cox model, adjustment for demogrpahic variables",
                             "Region-stratfied Cox model, adjustment for demogrpahic and clinical variables")))
  
  bind_rows(
    plot_data %>% mutate(outcome = "postest"),
    plot_data %>% mutate(outcome = "covidadmitted"),
    plot_data %>% mutate(outcome = "coviddeath"),
    plot_data %>% mutate(outcome = "death")
  ) %>%
    mutate(across(outcome, 
                  factor,
                  levels = c("postest", "covidadmitted", "coviddeath", "death"))) %>%
    ggplot(aes(x = comparison, y = estimate, colour = model)) +
    geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.1)) +
    geom_point(position = position_dodge(width = 0.1)) +
    geom_hline(aes(yintercept=1), colour='grey') +
    facet_wrap(~outcome) +
    # scale_y_log10(
    #   breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    #   limits = c(0.05, max(1, (plot_data$upper))),
    #   oob = scales::oob_keep,
    #   sec.axis = sec_axis(
    #     ~(1-.),
    #     name="Effectiveness",
    #     breaks = c(-4, -1, 0, 0.5, 0.80, 0.9, 0.95, 0.98, 0.99),
    #     labels = function(x) {formatpercent100(x, 1)}
    #   )
    # ) +
    scale_colour_brewer(name = NULL, type="qual", palette="Set2", guide=guide_legend(ncol=1))+
    labs(
      y="Hazard ratio",
      x="days since second dose",
      colour=NULL,
      title=glue("JCVI group {group}"),
      subtitle=glue("Two doses of {b} vs unvaccinated")
    ) +
    theme_bw()+
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text.x = element_text(size=8),
      
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = "bottom"
    ) 
  
}




  ##########

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('lubridate')
library('survival')
library('splines')
library('gtsummary')


## Import custom user functions from lib

# source(here("analysis", "lib", "utility_functions.R"))
# source(here("analysis", "lib", "redaction_functions.R"))
# source(here("analysis", "lib", "survival_functions.R"))


## import metadata ----
var_labels <- read_rds(here("output", "data", "metadata_labels.rds"))

list_formula <- read_rds(here("output", "data", "metadata_formulas.rds"))
list2env(list_formula, globalenv())

if(censor_seconddose==1){
  postvaxcuts<-postvaxcuts12
  lastfupday<-lastfupday12
} else {
  postvaxcuts <- postvaxcuts20
  if(outcome=="postest") postvaxcuts <- postvaxcuts20_postest
  lastfupday <- lastfupday20
}


# Import processed data ----

# import models ----

## report models ----

#data_cox <- read_rds(here("output", "models", outcome, timescale, "modelcox_data.rds"))

tidy_cox <- model_tidy
# <- read_rds(here("output", "models", outcome, timescale, censor_seconddose, "modelcox_tidy.rds"))

if(timescale == "calendar"){
  effectscox <- tidy_cox %>%
    filter(str_detect(term, fixed("unvax"))) %>%
    mutate(
      term=str_replace(term, pattern=fixed("vax1_az:timesincevax_pw"), ""),
      term=if_else(label=="vax1_az", paste0(postvaxcuts[1]+1,"-", postvaxcuts[2]), term),
      term=fct_inorder(term),
    )
}

if(timescale=="timesincevax"){
  effectscox <- tidy_cox %>%
    filter(str_detect(term, fixed("timesincevax_pw")) | str_detect(term, fixed("vax1_az"))) %>%
    mutate(
      term=str_replace(term, pattern=fixed("vax1_az:strata(timesincevax_pw)"), ""),
      term=fct_inorder(term),
    )
}

effectscox <- effectscox %>%
  mutate(
    term_left = as.numeric(str_extract(term, "^\\d+"))-1,
    term_right = as.numeric(str_extract(term, "\\d+$"))-1,
    term_right = if_else(is.na(term_right), lastfupday, term_right),
    term_midpoint = term_left + (term_right+1-term_left)/2
  )

write_csv(effectscox, path = here("output", "models", outcome, timescale, censor_seconddose, glue::glue("reportcox_effects.csv")))
write_rds(effectscox, path = here("output", "models", outcome, timescale, censor_seconddose, glue::glue("reportcox_effects.rds")))

plotcox <-
  ggplot(data = effectscox) +
  geom_point(aes(y=exp(estimate), x=term_midpoint, colour=model_name), position = position_dodge(width = 1.8))+
  geom_linerange(aes(ymin=exp(conf.low), ymax=exp(conf.high), x=term_midpoint, colour=model_name), position = position_dodge(width = 1.8))+
  geom_hline(aes(yintercept=1), colour='grey')+
  scale_y_log10(
    breaks=c(0.25, 0.33, 0.5, 0.67, 0.80, 1, 1.25, 1.5, 2, 3, 4),
    sec.axis = dup_axis(name="<--  favours Pfizer  /  favours AZ  -->", breaks = NULL)
  )+
  scale_x_continuous(breaks=unique(effectscox$term_left), limits=c(min(effectscox$term_left), max(effectscox$term_right)+1), expand = c(0, 0))+
  scale_colour_brewer(type="qual", palette="Set2", guide=guide_legend(ncol=1))+
  labs(
    y="Hazard ratio",
    x="Time since first dose",
    colour=NULL,
    title=glue::glue("AZ versus Pfizer, by time since first-dose")
  ) +
  theme_bw()+
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = "bottom"
  ) +
  NULL
plotcox
## save plot

write_rds(plotcox, path=here("output", "models", outcome, timescale, censor_seconddose, glue("reportcox_effectsplot.rds")))
ggsave(filename=here("output", "models", outcome, timescale, censor_seconddose, glue("reportcox_effectsplot.svg")), plotcox, width=20, height=15, units="cm")
ggsave(filename=here("output", "models", outcome, timescale, censor_seconddose, glue("reportcox_effectsplot.png")), plotcox, width=20, height=15, units="cm")
