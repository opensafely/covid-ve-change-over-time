################################################################################

# This script:
# - plots the vaccine effectiveness estimates for all outcomes and models for each comparison


################################################################################
library(tidyverse)
library(RColorBrewer)
library(glue)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  plot <- "BNT162b2andChAdOx" # "BNT162b2"  "ChAdOx" "BNT162b2andChAdOx" "BNT162b2vsChAdOx"
  
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

if (any(comparison %in% c("BNT162b2", "ChAdOx"))) {
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
                           "Stratfied Cox model, adjustment for demographic variables",
                           "Stratfied Cox model, adjustment for demographic and clinical variables"))) %>%
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


y_axis_label <- "Hazard ratio\n<--  favours vaccine  |  favours no vaccine  -->"
colour_var <- "model"
shape_var <- "subgroup"
caption_string <- "Stratification variables are: JCVI group, eligibility date for first dose of vaccination, geographical region."
shape_vals <- c("square", "circle", "triangle")[which(c("02", "03-10", "11-12") %in% subgroups)]
set2_vals <- brewer.pal(n = 4, name = "Set2")
colour_vals <- set2_vals[c(1,4)]
colour_name <- NULL
position_dodge_val <- 0.2
y_upper <- 2

if (plot %in% c("BNT162b2", "ChAdOx")) {
  title_string <- glue("Two doses of {plot} vs unvaccinated")
  if (plot %in% "BNT162b2") position_dodge_val <- 0.6
} else if (plot %in% "BNT162b2vsChAdOx") {
  title_string <- "Two doses of BNT162b2 vs two doses of ChAdOx"
  y_axis_label <- "Hazard ratio\n<--  favours BNT162b2  |  favours ChAdOx  -->"
  y_upper <- 10
} else if (plot %in% "BNT162b2andChAdOx") {
  colour_var <- "comparison"
  colour_vals <- set2_vals[c(2,3)]
  colour_name <- "Vaccine"
  title_string <- "Two doses of vaccine vs unvaccinated"
  caption_string <- str_wrap(
    "Hazard ratios estimated using a stratified cox model adjusted for demographic and clinical variables (stratification variables are: JCVI group, eligibility date for first dose of vaccination, geographical region)",
    150)
}




# set the colours corresponding to subgroups / brands
# shapes: triangle = unadjusted, circle = adjusted

plot_res <- plot_data %>%
  ggplot(aes(x = k, y = estimate, shape = !! sym(shape_var), colour = !! sym(colour_var))) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = position_dodge_val)) +
  geom_point(position = position_dodge(width = position_dodge_val)) +
  facet_wrap(~outcome, nrow=2, ncol=2)  +
  scale_y_log10(
    name = y_axis_label,
    breaks = c(0.00, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
    limits = c(0.01, max(1, y_upper)),
    oob = scales::oob_keep
  ) +
  # guides(shape = guide_legend(order = 1), 
  #        colour = guide_legend(order = 2)) +
  scale_shape_manual(name = "JCVI group(s)",  guide=guide_legend(nrow=1, order = 1), values = shape_vals) +
  scale_colour_manual(name = colour_name, guide=guide_legend(ncol=1, order = 2), values = colour_vals) +
  labs(
    x = "days since second dose",
    colour = NULL,
    title = title_string,
    caption = caption_string
  ) +
  theme_bw() +
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
       filename = here::here("output",  "images", glue("plot_res_{plot}.png")),
       width=20, height=15, units="cm")
