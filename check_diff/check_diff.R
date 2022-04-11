# check difference in estimates after corresction made to models

library(tidyverse)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

subgroup_width <- 25

outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds")
)
outcomes <- outcomes[outcomes != "covidemergency"]
old_names <- names(outcomes)
outcomes <- unname(outcomes)
new_names <- str_remove(old_names, "\\s\\(APCS\\)")
new_names <- str_remove(new_names, "SARS-CoV-2 ")
new_names <- str_replace(new_names, "COVID-19", "C-19")
new_names <- str_wrap(new_names, 15)

estimates_old <- readr::read_csv(here::here("release20220226", "estimates_all.csv")) %>%
  mutate(across(subgroup, factor, labels = subgroups[c(2,4,3,1)])) %>%
  mutate(across(subgroup, as.character)) %>%
  mutate(across(subgroup, factor, levels = subgroups, labels = str_wrap(subgroups, subgroup_width))) %>%
  filter(variable == "k", !reference_row) %>%
  mutate(across(model, ~if_else(.x == "unadjusted1", "unadjusted", "adjusted"))) %>%
  mutate(across(comparison, ~str_replace(.x, "ChAdOx", "ChAdOx1"))) %>%
  select(subgroup, comparison, outcome, model, label, estimate, conf.low, conf.high) %>%
  arrange(subgroup, comparison, outcome, model, label)
  

estimates_new <- readr::read_csv(here::here("release_20220401", "estimates_all.csv")) %>%
  mutate(across(subgroup, factor, levels = 1:4, labels = str_wrap(subgroups, subgroup_width))) %>%
  filter(variable == "k", !reference_row, outcome != "covidemergency") %>%
  select(subgroup, comparison, outcome, model, label, estimate, conf.low, conf.high) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) %>%
  arrange(subgroup, comparison, outcome, model, label)


estimates <- bind_rows(
  estimates_old %>%
    mutate(corrected = FALSE),
  estimates_new %>%
    mutate(corrected = TRUE)
) %>%
  mutate(across(model, factor, levels = c("unadjusted", "adjusted"))) %>%
  mutate(across(outcome, factor, levels = outcomes, labels = new_names))

gg_color_hue <- function(n, transparency = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
}

for (plot_comparison in c("BNT162b2", "ChAdOx1", "both")) {
  
  position_dodge_val <- 0.8
  
  if (plot_comparison == "BNT162b2") {
    col <- gg_color_hue(n=2)[1]
    plot_title <- str_c(plot_comparison, " vs unvaccinated")
  } else if (plot_comparison == "ChAdOx1") {
    col <- gg_color_hue(n=2)[2]
    plot_title <- str_c(plot_comparison, " vs unvaccinated")
  } else if (plot_comparison == "both") {
    col <- gg_color_hue(n=3)[2]
    plot_title <- "BNT162b2 vs ChAdOx1"
  }
  
  plot_out <- estimates %>%
    filter(comparison == plot_comparison) %>%
    ggplot(
      aes(x = label, y = estimate, alpha = model, shape = corrected)
    ) +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high), 
      position = position_dodge(width = position_dodge_val),
      colour = col
    ) +
    geom_point(
      position = position_dodge(width = position_dodge_val),
      colour = col
    ) +
    facet_grid(outcome ~ subgroup, scales = "free") +
    scale_alpha_manual(
      values = c("unadjusted" = 0.3, "adjusted" = 1)
    ) +
    # scale_color_manual(
    #   name = NULL,
    #   values = colour_pal
    # ) +
    scale_y_log10(
      name = "Hazard ratio",
      # breaks = c(0.00, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1),
      # limits = c(y_lower, y_upper),
      oob = scales::oob_keep
    ) +
    labs(
      x = "Weeks since second dose",
      title = plot_title
    ) +
    guides(alpha = guide_legend(order = 1), 
           shape = guide_legend(order = 2)) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=8),
      
      axis.title.x = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0, size = 12, face = "bold"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      plot.margin = margin(t=10, r=15, b=10, l=10),
      
      legend.position = "bottom"
      
    )
  
  ggsave(
    plot = plot_out,
    filename = here::here("release_20220401", glue("check_{plot_comparison}.png")),
    width = 24, height = 17, units = "cm"
  )
  
}


