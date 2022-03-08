library(tidyverse)
library(glue)

# path for released results
release_folder <- here::here("release20220226")

eligibility_count_ab <- readr::read_csv(file.path(release_folder, "eligibility_count_ab.csv"))
eligibility_count_ab$criteria <- NA_character_
eligibility_count_ab[1:9,]$criteria <- "a"
eligibility_count_ab[10:13,]$criteria <- "b"

eligibility_count_cde <- readr::read_csv(file.path(release_folder, "eligibility_count_cde.csv"))
eligibility_count_cde$criteria <- NA_character_
eligibility_count_cde[1,]$criteria <- "c"
eligibility_count_cde[2:4,]$criteria <- "e"
eligibility_count_cde[5,]$criteria <- "d"
eligibility_count_cde[6:8,]$criteria <- "e"

eligibility_count_p1 <- readr::read_csv(file.path(release_folder, "eligibility_count_p1.csv"))

eligibility_count_numeric <- bind_rows(
  eligibility_count_ab,
  eligibility_count_cde,
  eligibility_count_p1 %>% mutate(criteria = "p1")
) %>%
  mutate(arm = str_remove(str_extract(description, "^\\w+:"), ":"))

n_groups <- eligibility_count_numeric %>%
  mutate(across(criteria,
                ~if_else(
                  is.na(arm),
                  .x,
                  str_c(arm, .x, sep = "_")
                ))) %>%
  group_by(criteria) %>%
  summarise(
    min_n = min(n, na.rm=TRUE), 
    max_n = min(n, na.rm=TRUE), 
    .groups = "keep"
    ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = criteria,
    values_from = c(min_n, max_n)
  )

eligibility_count <- eligibility_count_numeric %>%
  mutate(across(starts_with("n"),
                ~ str_c("(n = ", scales::comma(.x, accuracy=1), ")")))

# create tibble for results
ncol <- 5
nrow <- 100
flow <- matrix("", nrow=nrow, ncol=ncol)

prior_covid <- str_c("(n = ", scales::comma(sum(eligibility_count_numeric[6:8,]$n_removed), accuracy = 1), ")")

# (not meeting inclusion B or D) - (meeting inclusion and exclusion A) - (meeting inclusion and exclusion B) - (meeting inclusion D)
not_b_or_d <- n_groups$min_n_a - n_groups$min_n_b - n_groups$min_n_vax_d 
not_b_or_d <- str_c("(n =", scales::comma(not_b_or_d, accuracy = 1), ")")

not_c <- eligibility_count_numeric[13,]$n - eligibility_count_numeric[18,]$n
not_c <- str_c("(n =", scales::comma(not_c, accuracy = 1), ")")

i <- 1
# A inc.
flow[i,3] <- glue("Meeting inclusion criteria A: assigned to JCVI groups 2-12 and registered for at least 1 year {eligibility_count[3,]$n}")
i <- i+2
# A ex.
flow[i,5] <- glue("Meeting exclusion criteria A: aged >120 {eligibility_count[4,]$n_removed}; missing  sex, ethnicity, region, or IMD {eligibility_count[5,]$n_removed}; evidence of prior COVID-19 {prior_covid}; in care home {eligibility_count[9,]$n_removed}")
i <- i+2
# A remaining
flow[i,3] <- glue("Remaining: {eligibility_count[9,]$n}")
i <- i+2
# Not B or D
flow[i,5] <- glue("Not meeting inclusion criteria B or D: {not_b_or_d}")
i <- i+2
# B inc.
flow[i,2] <- glue("Meeting inclusion criteria B: second dose of BNT162b2 or ChAdOx received between 6 and 14 weeks after eligibility date {eligibility_count[10,]$n}")
i <- i+2
# B ex
flow[i,1] <- glue("Meeting exclusion criteria B: first dose received before eligibility date {eligibility_count[11,]$n_removed}; <6 or >=14 weeks between first and second dose {eligibility_count[12,]$n_removed}; healthcare worker {eligibility_count[13,]$n_removed}")
i <- i+2
# Not C
flow[i,1] <- glue("Not meeting inclusion criteria C: {not_c}") 
i <- i+2
# C inc.
flow[i,2] <- glue("Meeting inclusion criteria C: first dose received during SVP {eligibility_count[18,]$n}")
# D inc
flow[i,4] <- glue("Meeting inclusion criteria D: unvaccinated at start of SVP {eligibility_count[14,]$n}")
i <- i+2
# E ex vax
flow[i,1] <- glue("Meeting exclusion citeria E: evidence of COVID-19 before start of SVP {eligibility_count[19,]$n_removed}; resident in care home before start of SVP {eligibility_count[20,]$n_removed}; end of life care before start of SVP {eligibility_count[21,]$n_removed}")
# E ex unvax
flow[i,5] <- glue("Meeting exclusion citeria E: evidence of COVID-19 before start of SVP {eligibility_count[15,]$n_removed}; resident in care home before start of SVP {eligibility_count[16,]$n_removed}; end of life care before start of SVP {eligibility_count[17,]$n_removed}")
i <- i+2
# P1 exclusions
flow[i,1] <- glue("Occurrence of the following before the start of period 1: death {eligibility_count[27,]$n_removed}; de-registration {eligibility_count[28,]$n_removed}; subsequent dose {eligibility_count[29,]$n_removed}")
flow[i,5] <- glue("Occurrence of the following before the start of period 1: death {eligibility_count[23,]$n_removed}; de-registration {eligibility_count[24,]$n_removed}; subsequent dose {eligibility_count[25,]$n_removed}")
# remaining
i <- i+2
flow[i,2] <- glue("Eligible for vaccinated group in period 1 {eligibility_count[29,]$n}")
flow[i,4] <- glue("Eligible for unvaccinated group in period 1 {eligibility_count[25,]$n}")

readr::write_delim(
  as.data.frame(flow),
  file = file.path(release_folder, "flow.csv"),
  delim = "\t"
)

