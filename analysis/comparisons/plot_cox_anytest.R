################################################################################

# This script:
# - plots the vaccine effectiveness estimates for all outcomes and models for each comparison


################################################################################
library(tidyverse)
library(RColorBrewer)
library(glue)
library(cowplot)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  plot <- "BNT162b2" # "BNT162b2"  "ChAdOx" "BNT162b2andChAdOx" "BNT162b2vsChAdOx"
  
} else {
  
  plot <- args[[1]]
  
}

################################################################################
fs::dir_create(here::here("output", "models_cox", "images"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "second_vax_period", "data", "second_vax_period_dates.rds")) 

outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)
plot_outcomes <- outcomes[outcomes=="anytest"]

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
# subgroups <- c(subgroups, "all")
subgroup_labels_full <- seq_along(subgroups)

# min and max follow-up dates per subgroup
min_and_max_fu_dates <- readr::read_rds(
  here::here("output", "lib", glue("min_and_max_fu_dates.rds")))

gg_color_hue <- function(n, transparency = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
}

################################################################################
arm1 <- if_else(plot =="ChAdOx", "ChAdOx", "BNT162b2")
arm2 <- if_else(plot == "BNT162b2vsChAdOx", "ChAdOx", "unvax")
arm3 <- if_else(plot == "BNT162b2andChAdOx", "ChAdOx", NA_character_)

data_tests <- readr::read_rds(
  here::here("output", "data", "data_tests.rds")) %>%
  select(patient_id, starts_with("pos_rate")) %>% 
  pivot_longer(
    cols = -patient_id,
    names_pattern = "pos_rate_(\\d)",
    names_to = "comparison",
    values_drop_na = TRUE
  ) 

data_comparisons <- local({
  
  data_arm1 <-  readr::read_rds(
    here::here("output", "comparisons", "data", glue("data_comparisons_{arm1}.rds"))) %>%
    select(patient_id, comparison, arm, subgroup)
  subgroups_1 <- unique(as.character(data_arm1$subgroup))
  
  data_arm2 <-  readr::read_rds(
    here::here("output", "comparisons", "data", glue("data_comparisons_{arm2}.rds"))) %>%
    select(patient_id, comparison, arm, subgroup)
  subgroups_2 <- unique(as.character(data_arm2$subgroup))
  subgroups <- intersect(subgroups_1, subgroups_2)
  
  data <- bind_rows(data_arm1, data_arm2)
  
  if (!is.na(arm3)) {
    data_arm3 <-  readr::read_rds(
      here::here("output", "comparisons", "data", glue("data_comparisons_{arm3}.rds"))) %>%
      select(patient_id, comparison, arm, subgroup)
    subgroups_3 <- unique(as.character(data_arm3$subgroup))
    subgroups <- intersect(subgroups, subgroups_3)
    data <- bind_rows(data, data_arm3)
  }
  
  data %>%
    filter(subgroup %in% subgroups) %>%
    mutate(across(comparison, as.character))
  
})

data_posrate <- data_comparisons %>%
  inner_join(data_tests, by = c("patient_id", "comparison"))

##############################################################################

if (plot %in% "BNT162b2") {
  subgroup_labels <- subgroup_labels_full
} else {
  subgroup_labels <- subgroup_labels_full[-which(subgroups == "18-39 years")]
}

if (plot == "BNT162b2andChAdOx") {
  comparisons <- c("BNT162b2", "ChAdOx")
} else if (plot == "BNT162b2vsChAdOx") {
  comparisons <- "both"
} else {
  comparisons <- plot
}

################################################################################
models <- c("0","2")

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

model_tidy_tibble <- bind_rows(
  model_tidy_list[sapply(model_tidy_list, function(x) is_tibble(x))]
) %>%
  mutate(across(outcome, factor, levels = plot_outcomes, labels = names(plot_outcomes)))

plot_fun <- function(
  plot_subgroup,
  plot_model,
  plot_comparison
) {
  
  args23 <- c("plot_model", "plot_comparison")
  
  # which one has length > 1
  lengths <- sapply(
    list(
      plot_model,
      plot_comparison
    ),
    function(x) length(x)
  )
  
  if (!(length(plot_subgroup) > 1 & (sum(lengths > 1) == 1)))
    stop("Plot_subgroup and one of plot_model and plot_comparison must have length>1")
  
  colour_var <- str_remove(args23[lengths > 1], "plot_")
  colour_var_length <- lengths[lengths > 1]
  
  # scale for x-axis
  K <- study_parameters$max_comparisons
  ends <- seq(2, (K+1)*4, 4)
  starts <- ends + 1
  weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")
  
  plot_data <- model_tidy_tibble %>% 
    filter(
      model %in% plot_model,
      subgroup %in% plot_subgroup,
      comparison %in% plot_comparison
    ) %>%
    filter(str_detect(term, "^comparison")) %>%
    mutate(
      k = as.integer(str_remove(str_extract(term, "comparison_\\d"), "comparison_"))
    )
  
  expanded_data <- tibble(
    outcome = character(),
    k = integer()
  )
  for (o in names(plot_outcomes)) {
    expanded_data <- expanded_data %>%
      bind_rows(tibble(
        outcome = rep(o, each = K),
        k = 1:K
      ))
  }
  
  # y-axis label
  y_axis_label <- "Hazard ratio\n<--  more testing in unvaccinated  |  more testing in vaccinated  -->"
  alpha_unadj <- 0.3
  # colour palette and name for colour legend
  if (colour_var == "model") {
    if (plot_comparison == "both") {
      
      palette_unadj <- gg_color_hue(3, transparency = alpha_unadj)
      palette_adj <- gg_color_hue(3, transparency = 1)
      i <- 2 # green
      
    } else {
      
      palette_unadj <- gg_color_hue(2, transparency = alpha_unadj)
      palette_adj <- gg_color_hue(2, transparency = 1)
      i <- case_when(
        plot_comparison %in% "BNT162b2" ~ 1,  # red 
        plot_comparison %in% "ChAdOx" ~ 2, # blue
        TRUE ~ NA_real_
      )
      
    }
    
    palette <- c(palette_unadj[i], palette_adj[i])
    colour_name <- NULL
  } else if (colour_var == "comparison") {
    # use ggplot palette so that this matches previous figures coloured by brand
    palette <- gg_color_hue(2)
    colour_name <- NULL
  }
  
  # spacing of points on plot
  position_dodge_val <- 0.6
  # upper limit for y-axis
  y_upper <- 10
  y_lower <- 0.1
  # plot caption
  caption_string <- if_else(
    colour_var == "model",
    "Stratification variables are: JCVI group, eligibility date for first dose of vaccination, geographical region.",
    "Hazard ratios estimated using a stratified Cox model adjusted for demographic and clinical variables (stratification variables are: JCVI group, eligibility date for first dose of vaccination, geographical region)"
  )
  x_title <- "weeks since second dose"
  
  subtitle_string <- str_c("Outcome: ", names(plot_outcomes),"\n ")
  
  if (plot %in% c("BNT162b2", "ChAdOx")) {
    
    title_string <- glue("Two doses of {plot} vs unvaccinated")
    
  } else if (plot %in% "BNT162b2vsChAdOx") {
    
    title_string <- "Two doses of BNT162b2 vs two doses of ChAdOx"
    y_axis_label <- "Hazard ratio\n<--  more testing in ChAdOx  |  more testing in BNT162b2  -->"
    
  } else if (plot %in% "BNT162b2andChAdOx") {
    
    title_string <- "Two doses of vaccine vs unvaccinated"
    
  }
  
    plot_height <- 16
    plot_width <- 15  
    legend_width <- 40
    caption_width <- 110
    theme_legend <- function(...) {
      theme(
        legend.position = "bottom",
        legend.key.size = unit(0.8, "cm"),
        ...
      )
    }
    ##
 
  plot_data2 <- expanded_data %>% 
    mutate(across(outcome, factor, levels = names(plot_outcomes))) %>%
    left_join(plot_data, by = c("outcome", "k")) %>%
    mutate(across(subgroup,
                  factor,
                  levels = subgroup_labels_full,
                  labels = sapply(subgroups, str_wrap, width=legend_width))) %>%
    left_join(min_and_max_fu_dates, by = "subgroup") %>%
    mutate(order = k) %>%
    mutate(across(k, 
                  factor, 
                  levels = 1:K,
                  labels = weeks_since_2nd_vax)) %>%
    group_by(subgroup) %>%
    mutate(
      min_k = min(order),
      max_k = max(order)
      ) %>%
    mutate(across(k,
                  ~ case_when(
                    order %in% min_k
                    ~ str_c(.x, "\n \nFollow-up from\n", min_fu),
                    order %in% max_k
                    ~ str_c(.x, "\n \nLatest follow-up\n", max_fu),
                    TRUE ~ as.character(.x)))) %>%
    mutate(across(model,
                  factor,
                  levels = models,
                  labels = sapply(c("Stratfied Cox model, no further adjustment", 
                                    "Stratfied Cox model, adjustment for demographic and clinical variables"),
                                  str_wrap, width=legend_width))) 
  
  plot_data2 %>% 
    filter(max_k < K) %>%
    distinct(subgroup, k) %>%
    mutate(k_new = )
  # pivot wider?
  
  
  plot_anytest <- plot_data2 %>%
    ggplot(aes(x = reorder(k,order), y = estimate, colour = !! sym(colour_var))) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = position_dodge_val)) +
    geom_point(position = position_dodge(width = position_dodge_val)) +
    facet_wrap(~subgroup, ncol=2, scales = "free_x")  +
    scale_y_log10(
      name = y_axis_label,
      breaks = c(0.00, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
      limits = c(y_lower, y_upper),
      oob = scales::oob_keep
    ) +
    scale_colour_manual(
      name = colour_name, 
      values = palette,
      na.translate = F) +
    labs(
      x = x_title,
      colour = NULL,
      title = title_string,
      subtitle = subtitle_string,
      caption = str_wrap(caption_string, caption_width)
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=8),
      
      axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 0, l = 0)),
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
      
    ) +
    theme_legend(
      legend.title = element_text(size = 10)
    )
  
  palette_posrate <- gg_color_hue(2, transparency = alpha_unadj)
  posrate_colours <- c(
    "BNT162b2" = palette_posrate[1], #red
    "ChAdOx" = palette_posrate[2], #blue
    "Unvaccinated" = "#C0C0C0" #grey
  )
  plot_data_posrate <- data_posrate %>%
    mutate(across(comparison, 
                  factor, 
                  levels = 1:K,
                  labels = weeks_since_2nd_vax)) %>%
    mutate(across(arm, ~if_else(.x%in%"unvax", "Unvaccinated", .x))) 
  
  plot_posrate <- plot_data_posrate %>%
    ggplot(aes(x = comparison, y = value, fill = arm)) +
    geom_boxplot(colour = "grey") +
    facet_wrap(~subgroup, nrow=2) +
    scale_fill_manual(
      name = NULL,
      values = posrate_colours[names(posrate_colours) %in% unique(plot_data_posrate$arm)]
      ) +
    labs(
      title = "SARS-CoV-2 test positivity rate",
      x = x_title,
      y = "crude rate*",
      caption = str_wrap("The crude SARS-CoV-2 test postivity rate was calculated as the total number of positive tests divided by the total number of tests in each comparison period. Individuals were removed from comparison periods during which they took zero tests.", caption_width)) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=8),
      
      axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 0, l = 0)),
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
      
    ) +
    theme_legend(
      legend.title = element_text(size = 10)
    )
  
  
  ic <- str_c(plot_subgroup, collapse = "")
  jc <- str_c(plot_model, collapse = "")
  
  # save the plots
  ggsave(plot = plot_anytest,
         filename = here::here("output", "models_cox", "images", glue("hr_anytest_{plot}_{ic}_{jc}.png")),
         width=plot_width, height=plot_height, units="cm")
  
  ggsave(plot = plot_posrate,
         filename = here::here("output", "models_cox", "images", glue("posrate_{plot}_{ic}_{jc}.png")),
         width=plot_width, height=plot_height, units="cm")
  
}

################################################################################

if (plot %in% c("BNT162b2", "ChAdOx", "BNT162b2vsChAdOx")) {
  j <- c("0", "2")
} else {
  j <- "2"
}

plot_fun(
  plot_subgroup = subgroup_labels,
  plot_model = j,
  plot_comparison = comparisons
)
  


