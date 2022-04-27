################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "BNT162b2"
  
} else{
  comparison <- args[[1]]
}

################################################################################
fs::dir_create(here::here("output", "models_cox", "images"))

################################################################################
# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)
if (comparison != "BNT162b2") {
  subgroup_labels <- subgroup_labels[subgroups != "18-39 years"]
}

################################################################################
for (i in subgroup_labels) {
  
  title_string <- glue("Subgroup: {subgroups[i]}")
  
  if (comparison %in% "BNT162b2" & i %in% 2) {
    plot_outcomes <- outcomes[outcomes != "coviddeath"]
  } else {
    plot_outcomes <- outcomes
  }
  
  
  # read summary data
  modelcox_summary <- lapply(
    unname(plot_outcomes),
    function(x)
      readr::read_rds(
        here::here("output", "models_cox", "data", glue("modelcox_tidy_{comparison}_{i}_{x}.rds"))
      ) %>% 
      mutate(outcome = x)
  )
  
  # create plot data
  plot_data <- bind_rows(
    modelcox_summary
  ) %>% 
    filter(variable != "k") %>%
    mutate(
      var_group = factor(case_when(
        str_detect(term, "^age") ~ "demographic",
        variable %in% "imd" ~ "demographic",
        variable %in% "sex" ~ "demographic",
        variable %in% "ethnicity" ~ "demographic",
        TRUE ~ "clinical"
      ),
      levels = c("demographic", "clinical"))) %>%
    mutate(across(outcome, factor, levels = unname(outcomes), labels = str_wrap(names(outcomes),18))) %>%
    arrange(var_group, term) %>%
    mutate(across(c(estimate, conf.low, conf.high), exp))
  
  # define order of variables for plot
  order <- plot_data %>%
    distinct(var_group, variable, term, label) %>%
    mutate(
      short_term = term,
      order = row_number()
    ) %>%
    mutate(across(short_term, ~str_remove(.x, " \\w+ deprived"))) %>%
    mutate(across(short_term, ~str_remove(.x, "ethnicity"))) %>%
    mutate(across(short_term, ~str_remove(.x, "bmi"))) %>%
    mutate(across(short_term, ~str_remove(.x, "TRUE"))) %>%
    mutate(across(short_term, ~str_trunc(.x, width = 15, side = "center"))) 
  
  # y_min <- min(plot_data$lower)
  # y_max <- max(plot_data$upper)
  
  # plot
  order %>%
    left_join(plot_data, by = c("term", "var_group")) %>%
    ggplot(aes(x = reorder(term,-order), y = estimate, colour = var_group)) +
    geom_hline(yintercept = 1, colour = "grey") +
    geom_point() +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
    facet_wrap(~ outcome, nrow=1, scales = "free_x") +
    scale_x_discrete(breaks = order$term, labels = order$short_term) +
    scale_y_log10(
      name = "hazard ratio",
      breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
      labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10"),
      # limits = c(min(0.1, y_min), max(10, y_max)),
      oob = scales::oob_keep
    ) +
    scale_color_discrete(guide = "none") +
    labs(
      subtitle = title_string,
      x = NULL
      ) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, size = 8),
      axis.text.y = element_text(size = 6)
    )
  ggsave(filename = here::here("output", "models_cox", "images", glue("coefs_{comparison}_{i}.png")),
         width=20, height=15, units="cm")
  
}
