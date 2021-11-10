### Function under development. May abandon so ignore for now.


# This is a function to create a frequency polygon that can be plotted by colour and/or faceted.
# The user can choose whether to output a plot or the data to generate a plot.
# A threshold can be set on the minimum number of samples in each bin.

redacted_stat_bin <- 
  function(
    .data, 
    xvar, # max length 1
    yvar_type=c("count","percent"),
    
    plot_type = c("freqpoly", "histogram"),
    
    facets = c("none", "grid", "wrap"),
    
    facet_vars = NULL, # max length 2, rows first then columns
    
    colour_var = NULL, # must be length 1
    
    binwidth=1, 
    
    plot = TRUE,
    
    threshold = 5, 
    
    ylim = c(0,NA),
    
    remove_tails=FALSE,
    
    ...) {
    
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(scales)
    library(tidyr)
    library(lubridate)
    
    # check_y <- deparse(substitute(yvar))
    # if (!(check_y %in% c('count', 'percent'))) stop('`yvar` must be `percent` or `count` (object not character)')
    # if (!(is.null(threshold) | is.numeric(threshold))) stop('`threshold` must be NULL or numeric.')
    # if (!is.numeric(binwidth)) stop('`binwidth` must be numeric.')
    # if (!(is.logical(facetwrap) && is.logical(facetgrid))) stop('`facetwrap` and `facetgrid` must both be logical.')
    
    # # include note at the bottom of the plot to state the binwidth
    # caption_string <- str_c('Binwidth = ', binwidth, if_else(is.null(units), '.', str_c(' ', units, '.')))
    caption_string <- "test"
    
    
    is_date <- unname(sapply(data_extract0 %>% select({{xvar}}), class))=="Date"
    
    data_in <- .data %>%
      filter(!is.na({{xvar}}))
    
    if (is_date) {
      data_in <- data_in %>%
        mutate(across({{xvar}}, as.integer))
    }
    
    rf <- 1/binwidth
    # lower and upper limits of data
    xvar_summary <- data_in %>% summarise(min = min({{xvar}}), 
                                          max = max({{xvar}}))
    xlower <- min(xvar_summary$min)
    xupper <- max(xvar_summary$max)
    
    # bins are closed on right, and I leave include.lowest(highest)=FALSE
    # in cut function, as I've +/- 0.5 on the max/min, so the lowest is not in the data
    # (this is consistent with stat_bin)
    # round lower (and upper) data to closest 0.5*binwidth
    xlower_b <- floor(xlower*rf)/rf - 0.5*binwidth
    xupper_b <- ceiling(xupper*rf)/rf + 0.5*binwidth
    
    my_seq <- signif(seq(xlower_b, xupper_b, binwidth),5)
    if (is_date) {
      # plot the points at start of the bins for date xvar
      labs_seq <- my_seq
    } else {
      # plot the points on the centre of the bins for numeric xvar
      labs_seq <- my_seq[-length(my_seq)] + 0.5*binwidth
    }
    
    plot_data <- data_in %>%
      as_tibble() %>%
      mutate(across(c({{colour_var}}, {{facet_vars}}), as.factor)) %>%
      mutate(across({{xvar}}, ~cut(.x,
                                   breaks = my_seq,
                                   # plot the points on the centre of the bins
                                   labels = labs_seq,
                                   right = TRUE))) %>%
      group_by(across(c({{xvar}}, {{facet_vars}}, {{colour_var}}))) %>%
      count() %>%
      ungroup() 
    
    # fill in zero counts
    expand_data <- plot_data %>%
      expand({{xvar}}, {{colour_var}}, {{facet_vars}})
    plot_data <- expand_data %>%
      left_join(plot_data) %>%
      mutate(across(n, ~if_else(is.na(.x), 0L, .x)))
    
    # remove zero counts if at tails of distribution
    if (remove_tails) {
      plot_data <- plot_data %>%
        group_by(across(c({{xvar}}, {{colour_var}}, {{facet_vars}}))) %>%
        arrange({{xvar}}, .by_group = TRUE) %>%
        mutate(sum_asc = cumsum(n)) %>%
        arrange(desc({{xvar}}), .by_group = TRUE) %>%
        mutate(sum_desc = cumsum(n)) %>%
        ungroup() %>%
        filter(sum_asc > 0, sum_desc > 0) %>%
        select(-sum_asc, - sum_desc) 
    }
    
    # if between zero and threshold, replace with threshold
    if (!is.null(threshold)) {
      plot_data <- plot_data %>%
        mutate(across(n, ~if_else((.x < threshold) & (.x > 0), 
                                  as.integer(threshold), 
                                  .x)))
    }
    
    # remove any duplicates
    plot_data <- plot_data %>%
      distinct()
    
    # convert xvar back to numeric or date
    if (is_date) {
      plot_data <- plot_data %>%
        mutate(across({{xvar}}, ~as.Date(.x, 
                                         format = "%Y-%m-%d", 
                                         origin = "1970-01-01"))) 
    } else {
      plot_data <- plot_data %>%
        mutate(across({{xvar}}, ~as.numeric(as.character(.x))))
    }
  
    if (plot) {
      # calculate percent and define yvar
      plot_data <- plot_data %>%
        group_by(across(c({{facet_vars}}, {{colour_var}}))) %>%
        # calculate percent of samples in each category
        mutate(percent = n/sum(n)) %>%
        ungroup() %>%
        rename("yvar" = yvar_type) %>%
        select({{xvar}}, {{facet_vars}}, {{colour_var}}, yvar)
      
      # create plots
      if (plot_type == "freqpoly") {
        plot <- plot_data %>%
          ggplot(aes(x = {{xvar}}, y = yvar, colour = {{colour_var}})) +
          geom_line() +
          lims(y = ylim) +
          labs(caption = caption_string)
      } else {
        plot <- plot_data %>%
          ggplot(aes(x = {{xvar}}, y = yvar, colour = {{colourvar}})) +
          geom_bar(stat = 'identity', width = binwidth) +
          lims(y = ylim) +
          labs(caption = caption_string)
      }
      
      # faceting
      if (facets == "grid") {
        plot <- plot +
          facet_grid(rows = vars({{facet_vars[1]}}),
                     cols = vars({{facet_vars[2]}}),
                     ...)
      } else if (facets == "wrap") {
        plot <- plot +
          facet_wrap(facets = vars({{facet_vars}}), ...)
      }
      
      # scales
      if (yvar_type == "percent") {
        plot <- plot +
          scale_y_continuous(labels = percent_format(), limits = ylim) +
          labs(y = NULL)
      } else {
        plot <- plot +
          lims(y = ylim) 
      }
      
      # theme
      # plot +
      #   theme_elsie()
      
      return(plot)
    } else {
      return(plot_data)
    }
  }
