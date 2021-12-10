

################################################################################
b <- "BNT162b2"

data_tte <- comparison_arms_data %>%
  filter(brand %in% b) %>%
  left_join(input_covs %>%
              select(patient_id, positive_test_date),
            by = "patient_id") %>%
  mutate(across(positive_test_date,
                ~ if_else(.x <= time_zero | .x > end_fu_date,
                          NA_Date_,
                          .x))) %>%
  group_by(brand, k) %>%
  mutate(origin = min(time_zero)) %>%
  ungroup() %>%
  mutate(start = as.integer(time_zero - origin),
         end = if_else(
           is.na(positive_test_date),
           as.integer(end_fu_date - origin),
           as.integer(positive_test_date - origin)),
         status = as.integer(!is.na(positive_test_date)))