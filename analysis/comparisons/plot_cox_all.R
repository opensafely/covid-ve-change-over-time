################################################################################
# This script:
# - combines estimates from all models into csv for release

################################################################################
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(glue)

################################################################################
fs::dir_create(here::here("output", "report", "data"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)
outcomes_order <- c(3,4,2,5)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)
subgroups_order <- c(4,1,3,2)

# define comparisons
comparisons <- c("BNT162b2", "ChAdOx", "both")

# min and max follow-up dates per subgroup
min_max_fu_dates <- readr::read_rds(
  here::here("output", "lib", glue("data_min_max_fu.rds"))) %>%
  mutate(across(ends_with("date"),
                ~ str_c(day(.x), " ", month(.x, label=TRUE))))

# scale for x-axis
K <- study_parameters$max_comparisons
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

################################################################################
model_tidy_list <- unlist(lapply(
  comparisons,
  function(x)
    unlist(lapply(
      subgroup_labels,
      function(y)
        lapply(
          unname(outcomes),
          function(z)
            try(
              readr::read_rds(
                here::here("output", "models_cox", "data", glue("modelcox_tidy_{x}_{y}_{z}.rds")
                )
              ) %>%
                mutate(comparison = x, subgroup = y, outcome = z)
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
) 

################################################################################

legend_width <- 15

plot_data <- model_tidy_tibble %>%
  filter(
    !reference_row,
    variable %in% "k",
    outcome != "anytest",
    model == 2
  ) %>%
  mutate(k=as.integer(label)) %>%
  group_by(subgroup, k) %>%
  mutate(k_missing = is.na(mean(estimate))) %>%
  ungroup() %>%
  filter(!(k_missing & subgroup==2)) %>%
  group_by(subgroup) %>%
  mutate(max_k = max(k)) %>%
  ungroup() %>%
  left_join(
    min_max_fu_dates %>%
      mutate(across(subgroup, ~subgroup_labels[subgroups == .x])), 
    by = "subgroup"
  ) %>%
  mutate(order1 = 10*k) %>%
  mutate(across(outcome,
                factor,
                levels = unname(outcomes[outcomes_order]),
                labels = str_wrap(names(outcomes[outcomes_order]), 10))) %>%
  mutate(k_labelled = k) %>%
  mutate(across(k_labelled, 
                factor, 
                levels = 1:K,
                labels = weeks_since_2nd_vax)) %>%
  mutate(across(model,
                factor,
                levels = 1:2,
                labels = sapply(c("Stratfied Cox model, no further adjustment", 
                                  "Stratfied Cox model, adjustment for demographic and clinical variables"),
                                str_wrap, width=legend_width))) %>%
  mutate(across(subgroup,
                factor,
                levels = subgroup_labels[subgroups_order],
                labels = str_wrap(subgroups[subgroups_order], 25))) %>%
  mutate(across(k_labelled,
                ~ case_when(
                  k == 1
                  ~ str_c(.x, "\n \nFollow-up\nfrom\n", min_fu_date),
                  k == max_k
                  ~ str_c(.x, "\n \nLatest\nfollow-up\n", max_fu_date),
                  TRUE ~ as.character(.x)))) %>%
  arrange(k) %>%
  group_by(k, subgroup, outcome, comparison) %>%
  mutate(order2 = row_number()) %>%
  ungroup() %>%
  mutate(order = order1 + order2)
    
# spacing of points on plot
position_dodge_val <- 0.6
# shape of points
point_shapes <- 21:24

# vaccine vs unvaccinated
plot_data %>%
  filter(comparison != "both") %>%
  ggplot(aes(
    x = reorder(k_labelled,order), 
    y = estimate, 
    colour = comparison, 
    fill = comparison, 
    shape = subgroup)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = position_dodge_val)) +
  geom_point(position = position_dodge(width = position_dodge_val)) +
  facet_grid(outcome ~ subgroup, switch = "y", scales = "free", space = "free_x") +
  scale_y_log10(
    name = "Hazard ratio\n<--   favours vaccine   |   favours no vaccine   -->",
    breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    limits = c(0.01, 5),
    oob = scales::oob_keep
  ) +
  labs(
    x = "weeks since second dose"
  ) +
  scale_colour_discrete(name = NULL) +
  scale_fill_discrete(guide = "none") +
  scale_shape_manual(guide = "none", values = point_shapes) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=8),
    
    axis.title.x = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = "bottom",
    legend.text = element_text(size=8)
  ) 

# brand comparison
plot_data %>%
  filter(comparison != "both") %>%
  ggplot(aes(
    x = reorder(k_labelled,order), 
    y = estimate, 
    colour = comparison, 
    fill = comparison, 
    shape = subgroup)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = position_dodge_val)) +
  geom_point(position = position_dodge(width = position_dodge_val)) +
  facet_grid(outcome ~ subgroup, switch = "y", scales = "free", space = "free_x") +
  scale_y_log10(
    name = "Hazard ratio\n<--   favours vaccine   |   favours no vaccine   -->",
    breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    limits = c(0.01, 5),
    oob = scales::oob_keep
  ) +
  labs(
    x = "weeks since second dose"
  ) +
  scale_colour_discrete(name = NULL) +
  scale_fill_discrete(guide = "none") +
  scale_shape_manual(guide = "none", values = point_shapes) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=8),
    
    axis.title.x = element_text(size = 8, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = "bottom",
    legend.text = element_text(size=8)
  ) 

# save the plot
ggsave(filename = here::here("output", "models_cox", "images", glue("hr_all.png")),
       width=20, height=15, units="cm")


# %>%
#   mutate(panel_order = case_when(
#     subgroup %in% subgroup_labels[4] & outcome %in% unname(outcomes[4]) ~ 1,
#     subgroup %in% subgroup_labels[4] & outcome %in% unname(outcomes[5]) ~ 2,
#     subgroup %in% subgroup_labels[1] & outcome %in% unname(outcomes[4]) ~ 3,
#     subgroup %in% subgroup_labels[1] & outcome %in% unname(outcomes[5]) ~ 4,
#     subgroup %in% subgroup_labels[4] & outcome %in% unname(outcomes[2]) ~ 5,
#     subgroup %in% subgroup_labels[4] & outcome %in% unname(outcomes[3]) ~ 6,
#     subgroup %in% subgroup_labels[1] & outcome %in% unname(outcomes[2]) ~ 7,
#     subgroup %in% subgroup_labels[1] & outcome %in% unname(outcomes[3]) ~ 8,
#     subgroup %in% subgroup_labels[3] & outcome %in% unname(outcomes[4]) ~ 9,
#     subgroup %in% subgroup_labels[3] & outcome %in% unname(outcomes[5]) ~ 10,
#     subgroup %in% subgroup_labels[2] & outcome %in% unname(outcomes[4]) ~ 11,
#     subgroup %in% subgroup_labels[2] & outcome %in% unname(outcomes[5]) ~ 12,
#     subgroup %in% subgroup_labels[3] & outcome %in% unname(outcomes[2]) ~ 13,
#     subgroup %in% subgroup_labels[3] & outcome %in% unname(outcomes[3]) ~ 14,
#     subgroup %in% subgroup_labels[2] & outcome %in% unname(outcomes[2]) ~ 15,
#     subgroup %in% subgroup_labels[2] & outcome %in% unname(outcomes[3]) ~ 16,
#     TRUE ~ NA_real_
#   )) %>%
#   mutate(across(subgroup, 
#                 factor, 
#                 levels = subgroup_labels[subgroups_order],
#                 labels = subgroups[subgroups_order])) %>%
#   mutate(across(outcome, 
#                 factor, 
#                 levels = unname(outcomes[outcomes_order]), 
#                 labels = names(outcomes[outcomes_order]))) %>%
#   mutate(panel = str_wrap(str_c(subgroup, outcome, sep = ", "),30)) 
# 
# panel_labels <- plot_data %>%
#   distinct(panel, panel_order) %>%
#   arrange(panel_order) %>%
#   select(panel) %>% unname() %>% unlist()
# 
# plot_data %>%
#   mutate(across(panel_order, levels = 1:16, labels = panel_labels)) %>%
#   ggplot(aes(x = label, y = estimate, colour = comparison)) +
#   geom_point() +
#   facet_wrap(~panel_order)
