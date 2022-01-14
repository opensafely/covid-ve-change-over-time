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

fs::dir_create(here::here("output", "models_cox", "images"))

# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "second_vax_period", "data", "second_vax_period_dates.rds")) 

outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroups <- c(subgroups, "all")
subgroup_labels <- seq_along(subgroups)

################################################################################

if (plot != "BNT162b2") {
  subgroup_labels <- subgroup_labels[-which(subgroups == "18-39")]
  subgroups <- subgroups[-which(subgroups == "18-39")]
} 

if (plot == "BNT162b2andChAdOx") {
  comparisons <- c("BNT162b2", "ChAdOx")
} else if (plot == "BNT162b2vsChAdOx") {
  comparison <- "both"
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
          outcomes,
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
) 




plot_fun <- function(
  plot_subgroup,
  plot_model,
  plot_comparison
) {
  
  args <- names(formals())
  
  # which one has length > 1
  lengths <- sapply(
    list(
      plot_subgroup,
      plot_model,
      plot_comparison
    ),
    function(x) length(x)
  )
  
  if (sum(lengths > 1) != 1)
    stop("Exactly one of plot_subgroup, plot_model and plot_comparison must have length>1")
  
  colour_var <- str_remove(args[lengths > 1], "plot_")
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
  for (o in outcomes) {
      expanded_data <- expanded_data %>%
        bind_rows(tibble(
          outcome = rep(o, each = K),
          k = 1:K
        ))
  }
  
  # y-axis label
  y_axis_label <- "Hazard ratio\n<--  favours vaccine  |  favours no vaccine  -->"
  # colour palette and name for colour legend
  if (colour_var == "subgroup") {
    set2vals <- brewer.pal(n = max(subgroup_labels), name = "Set2")[subgroup_labels]
    colour_name <- "Age range"
  } else {
    set2vals <- brewer.pal(n =  max(colour_var_length,3), name = "Set2")[1:colour_var_length]
    colour_name <- NULL
  }
  # spacing of points on plot
  position_dodge_val <- 0.6
  # upper limit for y-axis
  y_upper <- 2
  # plot caption
  caption_string <- if_else(
    colour_var == "model",
    "Stratification variables are: JCVI group, eligibility date for first dose of vaccination, geographical region.",
    "Hazard ratios estimated using a stratified cox model adjusted for demographic and clinical variables (stratification variables are: JCVI group, eligibility date for first dose of vaccination, geographical region)"
  )
  
  if (colour_var == "subgroup") {
    subtitle_string <- ""
  } else {
    subtitle_string <- str_c("(age range in plot: ", subgroups[plot_subgroup], ")")
  }
  
  if (plot %in% c("BNT162b2", "ChAdOx")) {
    
    title_string <- glue("Two doses of {plot} vs unvaccinated")
    
  } else if (plot %in% "BNT162b2vsChAdOx") {
    
    title_string <- "Two doses of BNT162b2 vs two doses of ChAdOx"
    y_axis_label <- "Hazard ratio\n<--  favours BNT162b2  |  favours ChAdOx  -->"
    y_upper <- 10
    
  } else if (plot %in% "BNT162b2andChAdOx") {
    
    title_string <- "Two doses of vaccine vs unvaccinated"
    
  }
  
  legend_width <- 40
  
  plot_res <- expanded_data %>% 
    mutate(across(outcome, factor)) %>%
    left_join(plot_data, by = c("outcome", "k")) %>%
    mutate(across(model,
                  factor,
                  levels = models,
                  labels = sapply(c("Stratfied Cox model, no further adjustment", 
                                    "Stratfied Cox model, adjustment for demographic and clinical variables"),
                                  str_wrap, width=legend_width))) %>%
    mutate(across(outcome,
                  factor,
                  levels = outcomes,
                  labels = c("Positive COVID-19 test",
                             "COVID-19 hospital admission",
                             "COVID-19 death",
                             "non-COVID-19 death",
                             "Any death"))) %>%
    mutate(across(subgroup,
                  factor,
                  levels = subgroup_labels,
                  labels = sapply(subgroups, str_wrap, width=legend_width))) %>%
    mutate(across(k, 
                  factor, 
                  levels = 1:K,
                  labels = weeks_since_2nd_vax)) %>%
    ggplot(aes(x = k, y = estimate, colour = !! sym(colour_var))) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = position_dodge_val)) +
    geom_point(position = position_dodge(width = position_dodge_val)) +
    facet_wrap(~outcome, nrow=3, ncol=2)  +
    scale_y_log10(
      name = y_axis_label,
      breaks = c(0.00, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
      limits = c(0.01, max(1, y_upper)),
      oob = scales::oob_keep
    ) +
    scale_colour_manual(
      name = colour_name, 
      guide=guide_legend(ncol=1), 
      values = set2vals,
      na.translate = F) +
    labs(
      x = "weeks since second dose",
      colour = NULL,
      title = title_string,
      subtitle = subtitle_string,
      caption = str_wrap(caption_string,110)
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
      
      legend.title = element_text(size = 10),
      legend.position = c(0.75, 0.15), # c(0,0) bottom left, c(1,1) top-right.
      legend.key.size = unit(0.8, "cm"),
      # legend.spacing.y = unit(2, 'cm'),
      legend.box="vertical"
    ) 
  
  return(plot_res)
  
}

if (plot %in% c("BNT162b2", "ChAdOx", "BNT162b2andChAdOx")) {
  
  plot_subgroups <- as.list(subgroup_labels)
  
  if (plot %in% c("BNT162b2", "ChAdOx")) {
    plot_subgroups <- splice(
      plot_subgroups,c(1:4)
    )
  }
  
  for (i in seq_along(plot_subgroups)) {
    
    if (length(plot_subgroups[[i]])==1 && plot %in% c("BNT162b2", "ChAdOx")) {
      j <- c("0", "2")
    } else {
      j <- "2"
    }
    
    plot_res <- plot_fun(
      plot_subgroup = plot_subgroups[[i]],
      plot_model = j,
      plot_comparison = comparisons
    )
    
    ic <- str_c(plot_subgroups[[i]], collapse = "")
    jc <- str_c(j, collapse = "")
    
    # save the plot
    ggsave(plot = plot_res,
           filename = here::here("output", "models_cox", "images", glue("plot_res_{plot}_{ic}_{jc}.png")),
           width=15, height=18, units="cm")
    
  }
  
}
  
  