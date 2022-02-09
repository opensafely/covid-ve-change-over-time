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
  plot <- "BNT162b2" # "BNT162b2"  "ChAdOx" "BNT162b2andChAdOx" "BNT162b2vsChAdOx"
  
} else {
  
  plot <- args[[1]]
  
}


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

gg_color_hue <- function(n, transparency = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
}

################################################################################

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
 
  
  
  
  plot_res <- expanded_data %>% 
    mutate(across(outcome, factor, levels = names(plot_outcomes))) %>%
    left_join(plot_data, by = c("outcome", "k")) %>%
    mutate(across(model,
                  factor,
                  levels = models,
                  labels = sapply(c("Stratfied Cox model, no further adjustment", 
                                    "Stratfied Cox model, adjustment for demographic and clinical variables"),
                                  str_wrap, width=legend_width))) %>%
    mutate(across(subgroup,
                  factor,
                  levels = subgroup_labels_full,
                  labels = sapply(subgroups, str_wrap, width=legend_width))) %>%
    mutate(across(k, 
                  factor, 
                  levels = 1:K,
                  labels = weeks_since_2nd_vax)) %>%
    ggplot(aes(x = k, y = estimate, colour = !! sym(colour_var))) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = position_dodge_val)) +
    geom_point(position = position_dodge(width = position_dodge_val)) +
    facet_wrap(~subgroup, ncol=2)  +
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
      x = "weeks since second dose",
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
  
  ic <- str_c(plot_subgroup, collapse = "")
  jc <- str_c(plot_model, collapse = "")
  
  # save the plot
  ggsave(plot = plot_res,
         filename = here::here("output", "models_cox", "images", glue("hr_anytest_{plot}_{ic}_{jc}.png")),
         width=plot_width, height=plot_height, units="cm")
  
  return(plot_res)
  
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
  


