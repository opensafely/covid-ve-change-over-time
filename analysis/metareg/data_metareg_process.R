library(tidyverse)

# read data
estimates_all <- readr::read_csv(
  here::here("output", "models_cox", "data", "estimates_all.csv"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))

# output folder
fs::dir_create(here::here("output", "metareg", "data"))

# prcoess
data_metareg <- estimates_all %>%
  filter(
    !reference_row,
    variable %in% "k",
    model %in% "unadjusted2"
    ) %>%
  mutate(model = "adjusted") %>%
  select(subgroup, comparison, outcome, model, k = label, estimate, conf.low, conf.high) %>%
  mutate(across(c(estimate, conf.low, conf.high), 
                ~ case_when(
                  is.na(.x) 
                  ~ mean(.x, na.rm=TRUE), 
                  (outcome %in% "coviddeath") & 
                    (subgroup %in% 1) &
                    ((comparison %in% c("BNT162b2", "both") & k < 4) |
                       comparison %in% "ChAdOx" & k < 2)
                  ~ NA_real_,
                  (comparison %in% c("BNT162b2", "both")) & 
                    (subgroup %in% 2) &
                    ((outcome %in% "covidadmitted" & k == 6) |
                       (outcome %in% "coviddeath") |
                       (outcome %in% "noncoviddeath" & k %in% c(2,5,6)))
                  ~ NA_real_,
                  (comparison %in% c("BNT162b2", "both")) & 
                    (subgroup %in% 3) &
                    ((outcome %in% "covidadmitted" & k != 4) |
                       (outcome %in% "coviddeath") |
                       (outcome %in% "noncoviddeath" & k %in% 6))
                  ~ NA_real_,
                  (comparison %in% c("ChAdOx")) & 
                    (subgroup %in% 3) &
                    ((outcome %in% "coviddeath" & k < 3))
                  ~ NA_real_,
                  (comparison %in% c("ChAdOx")) & 
                    (subgroup %in% 4) &
                    ((outcome %in% "coviddeath" & k == 1))
                  ~ NA_real_,
                  (comparison %in% c("BNT162b2")) & 
                    (subgroup %in% 4) &
                    ((outcome %in% "coviddeath" & k == 2))
                  ~ NA_real_,
                  (comparison %in% c("both")) & 
                    (subgroup %in% 4) &
                    ((outcome %in% "coviddeath" & k < 3))
                  ~ NA_real_,
                  TRUE 
                  ~ .x))) %>%
  mutate(across(subgroup, factor, levels = 1:4, labels = subgroups))

# save data
readr::write_csv(
  data_metareg,
  here::here("output", "metareg", "data", "data_metareg.csv"))
