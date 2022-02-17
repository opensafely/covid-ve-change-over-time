library(meta)


data_list <- model_tidy_tibble_arms %>%
  mutate(k = as.integer(str_extract(term, "\\d"))) %>%
  group_split(subgroup, comparison, outcome, model)

# assume linear relationship (period as continuous, intercept at period=1)
data_lin <- data_list[[1]] %>%
  mutate(k=k-1)

?metabin

# allow for linear relationship (period as factor, reference at period=1)