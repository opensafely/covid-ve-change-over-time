################################################################################

# This script:
# - plots the vaccine effectiveness estimates for all outcomes and models for each comparison


################################################################################
library(tidyverse)
library(glue)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  plot <- "BNT162b2andChAdOx" #"ChAdOx" "BNT162b2andChAdOx" "BNT162b2vsChAdOx"
  
} else {
   
  plot <- args[[1]]
  
}


################################################################################

fs::dir_create(here::here("output", "images"))

# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "lib", "second_vax_period_dates.rds")) 

outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)

################################################################################

subgroups <- "03-10"
if (plot == "BNT162b2") {
  subgroups <- c("02", "03-10", "11-12")
  comparison <- "BNT162b2"
}else if (plot == "ChAdOx") {
  comparison <- "ChAdOx"
} else if (plot == "BNT162b2andChAdOx") {
  comparison <- c("BNT162b2", "ChAdOx")
} else if (plot == "BNT162b2vsChAdOx") {
  comparison <- "both"
}

################################################################################

formatpercent100 <- function(x,accuracy){
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}

################################################################################

if (comparison %in% c("BNT162b2", "ChAdOx")) {
  plot_title <- glue("Two doses of {comparison} vs unvaccinated")
} else if (comparison %in% "both") {
  plot_title <- "Two doses of BNT162b2 vs two doses of ChAdOx"
}


models <- c("0","2")

model_tidy_list <- unlist(lapply(
  subgroups,
  function(x)
    unlist(lapply(
      comparison,
      function(y)
        lapply(
          outcomes,
          function(z)
            try(
              readr::read_rds(
                here::here("output", "models", glue("modelcox_summary_{x}_{y}_{z}.rds")
                )
              ) %>%
                mutate(subgroup = x, comparison = y)
            )
        )
    ),
    recursive = FALSE
    )
),
recursive = FALSE
)

model_tidy_tibble <- bind_rows(
  model_tidy_list[sapply(model_tidy_list, function(x) is_tibble(x))]
) %>%
  filter(str_detect(term, "^comparison")) 

if (plot == "BNT162b2andChAdOx") {
  model_tidy_tibble <- model_tidy_tibble %>%
    filter(model %in% "2")
}

K <- study_parameters$max_comparisons

ends <- seq(14, (K+1)*28, 28)
starts <- ends + 1
days_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

plot_data <- model_tidy_tibble %>% 
  mutate(
    k = factor(as.integer(str_remove(str_extract(term, "comparison_\\d"), "comparison_")),
                        labels = days_since_2nd_vax)
  ) %>%
  mutate(across(model,
                factor,
                levels = as.character(0:2),
                labels = c("Stratfied Cox model, no further adjustment", 
                           "Stratfied Cox model, adjustment for demogrpahic variables",
                           "Stratfied Cox model, adjustment for demogrpahic and clinical variables"))) %>%
  mutate(across(outcome,
                factor,
                levels = outcomes,
                labels = c("Positive COVID-19 test",
                           "COVID-19 hospital admission",
                           "COVID-19 death",
                           "Any death"))) %>%
  mutate(across(subgroup,
                factor,
                levels = c("02", "03-10", "11-12"),
                labels = c("2", "3-10", "11-12")))


if (plot %in% c("BNT162b2", "ChAdOx")) {
  colour_var <- "subgroup"
  title_string <- ""
} else if (plot %in% "BNT162b2andChAdOx") {
  colour_var <- "comparison"
  title_string <- ""
} else if (plot %in% "BNT162b2vsChAdOx") {
  colour_var <- "comparison"
  title_string <- ""
}



# set the colours corresponding to subgroups / brands
# shapes: triangle = unadjusted, circle = adjusted

plot_res <- plot_data %>%
  ggplot(aes(x = k, y = estimate, colour = !! sym(colour_var), shape = model)) +
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  facet_wrap(~outcome, nrow=2, ncol=2) 

if (plot %in% c("BNT162b2", "ChAdOx", "BNT162b2andChAdOx")) {
  plot_res <- plot_res +
    scale_y_log10(
      breaks = c(0.00, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
      limits = c(0.01, max(1, (plot_data$upper))),
      oob = scales::oob_keep
    )
} else if (plot %in% "BNT162b2vsChAdOx") {
  plot_res <- plot_res +
    scale_y_log10(
      breaks=c(0.25, 0.33, 0.5, 0.67, 0.80, 1, 1.25, 1.5, 2, 3, 4),
      sec.axis = dup_axis(name="HR in favour of\n<- ChAdOx / BNT162b2 ->", breaks = NULL)
    ) 
}

plot_res <- plot_res +
  scale_shape_discrete(name = NULL,  guide=guide_legend(ncol=1)) +
  scale_colour_brewer(name = "JCVI group(s)", type="qual", palette="Set2", guide=guide_legend(nrow=1)) +
  labs(
    y="Hazard Ratio (HR)",
    x="days since second dose",
    colour=NULL,
    title=plot_title,
    caption="Stratification variables are: JCVI group, eligibility date for first dose of vaccination, geographical region."
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
    
    legend.position = "bottom",
    legend.box="vertical"
  ) 

# save the plot
ggsave(plot = plot_res,
       filename = here::here("output",  "images", glue("plot_res_{comparison}.png")),
       width=20, height=15, units="cm")




