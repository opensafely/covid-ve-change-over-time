################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  group <- "02"
  
} else{
  group <- args[[1]]
}

################################################################################

fs::dir_create(here::here("output", glue("jcvi_group_{group}"), "images"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds")) %>%
  filter(jcvi_group %in% group, include) %>%
  distinct(brand, n_comparisons)

outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)


formatpercent100 <- function(x,accuracy){
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}

for (b in as.character(unique(second_vax_period_dates$brand))) {
  
  models <- as.character(0:2)

  model_tidy_list <- lapply(
    outcomes,
    function(x)
    try(
      readr::read_rds(
        here::here("output", glue("jcvi_group_{group}"), "models", glue("{b}_{x}_modelcox_summary.rds")
                   )
        )
      )
    )

  model_tidy_tibble <- bind_rows(
    model_tidy_list[sapply(model_tidy_list, function(x) is_tibble(x))]
    ) %>%
    filter(str_detect(term, "^armvax")) 
  
  
  K <- min(second_vax_period_dates$n_comparisons[second_vax_period_dates$brand == b])
  
  ends <- seq(14, (K+1)*28, 28)
  starts <- ends + 1
  days_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")
  
  plot_data <- model_tidy_tibble %>% 
    mutate(
      comparison = factor(as.integer(str_extract(term, "\\d")),
                          labels = days_since_2nd_vax)
    ) %>%
    mutate(across(model,
                  factor,
                  levels = as.character(0:2),
                  labels = c("Region-stratfied Cox model, no further adjustment", 
                             "Region-stratfied Cox model, adjustment for demogrpahic variables",
                             "Region-stratfied Cox model, adjustment for demogrpahic and clinical variables"))) %>%
    mutate(across(outcome,
                  factor,
                  levels = outcomes,
                  labels = c("Positive COVID-19 test",
                             "COVID-19 hospital admission",
                             "COVID-19 death",
                             "Any death")))
  
  plot_res <- plot_data %>%
    ggplot(aes(x = comparison, y = estimate, colour = model)) +
    geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.25)) +
    geom_point(position = position_dodge(width = 0.25)) +
    geom_hline(aes(yintercept=1), colour='grey') +
    facet_wrap(~outcome) +
    scale_y_log10(
      breaks = c(0.00, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
      limits = c(0.01, max(1, (plot_data$upper))),
      oob = scales::oob_keep,
      sec.axis = sec_axis(
        ~(1-.),
        name="Effectiveness (1 - HR)",
        breaks = c(-4, -1, 0, 0.5, 0.80, 0.9, 0.95, 0.98, 0.99, 1.00),
        labels = function(x) {formatpercent100(x, 1)}
      )
    ) +
    scale_colour_brewer(name = NULL, type="qual", palette="Set2", guide=guide_legend(ncol=1))+
    labs(
      y="Hazard Ratio (HR)",
      x="days since second dose",
      colour=NULL,
      title=glue("JCVI group {group}"),
      subtitle=glue("Two doses of {b} vs unvaccinated")
    ) +
    theme_bw()+
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text.x = element_text(size=6),
      
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
  
  # save the plot
  ggsave(plot = plot_res,
         filename = here::here("output", glue("jcvi_group_{group}"), "images", glue("plot_res_{b}.png")),
         width=20, height=15, units="cm")
  
}



