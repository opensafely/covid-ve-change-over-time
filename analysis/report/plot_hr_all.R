################################################################################

# This script:
# - plots the vaccine effectiveness estimates for all outcomes and models for each comparison


################################################################################
library(tidyverse)
library(RColorBrewer)
library(glue)

################################################################################
fs::dir_create(here::here("output", "report", "images"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)
plot_outcomes <- outcomes

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# define comparisons
comparisons <- c("BNT162b2", "ChAdOx", "both")

# define function for gg pallaette
gg_color_hue <- function(n, transparency = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
}

################################################################################
cat("read model data\n")
model_tidy_list <- unlist(lapply(
  comparisons,
  function(x)
    unlist(lapply(
      subgroup_labels,
      function(y)
        lapply(
          unname(plot_outcomes),
          function(z)
            try(
              readr::read_rds(
                here::here("output", "models_cox", "data", glue("modelcox_summary_{x}_{y}_{z}.rds")
                )
              ) %>%
                mutate(comparison = x, subgroup = y)
            )
        )
    ),
    recursive = FALSE
    )
),
recursive = FALSE
)

cat("bind model data\n")
model_tidy_tibble <- bind_rows(
  model_tidy_list[sapply(model_tidy_list, function(x) is_tibble(x))]
) %>%
  mutate(across(outcome, factor, levels = plot_outcomes, labels = names(plot_outcomes)))

cat("specify objects for plot\n")
# specify palette for plot
plot_pal <- gg_color_hue(n=3, transparency = 0.7)
# specify models to plot
models <- "2"
# specify scale for x-axis
K <- study_parameters$max_comparisons
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")
# specify dodge value
position_dodge_val <- 0.8
# specify orders in plot
subgroups_order <- c(4,1,3,2)
outcomes_order <- c(2,3,4,1,5)

cat("derive plot_data")
plot_data <- model_tidy_tibble %>%
  filter(
    str_detect(term, "^comparison"),
    model=="2"
    ) %>%
  mutate(
    k = str_extract(term, "\\d")
    ) %>%
  mutate(across(comparison, 
                ~if_else(.x %in% "both",
                         "BNT162b2 vs ChAdOx", 
                         str_c(.x, " vs unvaccinated")))) %>%
  mutate(across(subgroup, 
                factor, 
                levels = subgroup_labels[subgroups_order],
                labels = subgroups[subgroups_order])) %>%
  mutate(across(outcome,
                factor,
                levels = levels(model_tidy_tibble[["outcome"]])[outcomes_order])) %>%
  mutate(across(k, 
                factor, 
                levels = 1:K,
                labels = weeks_since_2nd_vax)) 

cat("create plot\n")
plot_res <- plot_data %>%
  ggplot(aes(x = k, y = estimate, colour = comparison, shape = subgroup)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_point(position = position_dodge(width = position_dodge_val)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 position = position_dodge(width = position_dodge_val)) +
  facet_wrap(~outcome, scales = "free_x", nrow=3) +
  scale_x_discrete(
    name = "weeks since second dose"
  ) +
  scale_y_log10(
    name = "<-- favours LHS -------- hazard ratio -------- favours RHS -->",
    breaks = c( 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    limits = c(0.01, 5),
    oob = scales::oob_keep#,
    # sec.axis = sec_axis(~(1-.), 
    #                     name="Effectiveness", 
    #                     breaks = c(-4, -1, 0, 0.5, 0.80, 0.9, 0.95, 0.98, 0.99), 
    #                     labels = scales::label_percent(1))
  ) +
  scale_colour_manual(
    name = "Comparison (LHS vs RHS)",
    values = c("BNT162b2 vs ChAdOx" = plot_pal[2],
               "BNT162b2 vs unvaccinated" = plot_pal[1],
               "ChAdOx vs unvaccinated" = plot_pal[3])
  ) +
  scale_shape_discrete(
    name = "Subgroup"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=8),
    
    axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    
    panel.spacing = unit(0.5, "lines"),
    
    plot.title = element_text(hjust = 0, size = 11),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = c(0.78, 0.15),
    legend.key.size = unit(0.5, "cm")
    
  ) 

cat("save plot\n")
ggsave(plot_res,
       filename = here::here("output", "report", "images", glue("plot_hr_all.png")),
       width=16, height=20, units="cm")
 

