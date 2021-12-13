

test <- data_comparison_arms %>%
  filter(comparison==1, brand == "BNT162b2") %>%
  left_join(input_covs %>%
              select(patient_id, positive_test_0_date),
            by = "patient_id") %>%
  mutate(origin = min(time_zero_date)) %>%
  mutate(across(c(time_zero_date, end_fu_date, positive_test_0_date),
                ~ as.integer(.x - origin)))
  
library(survival)
test_tmerge <- as_tibble(
  tmerge(
    data1 = test %>% select(patient_id),
    data2 = test,
    id = patient_id,
    tstart = time_zero_date,
    tstop = end_fu_date,
    
    postest = event(positive_test_0_date)
  )
)
  
as_tibble(pbc)
as_tibble(pbcseq)

temp <- as_tibble(subset(pbc, id <= 10, select=c(id:sex, stage))) # baseline data
pbc2 <- as_tibble(tmerge(temp, temp, id=id, endpt = event(time, status)))
pbc2 <- as_tibble(
  tmerge(pbc2, 
         subset(pbcseq, id <= 10), 
         id=id, 
         ascites = tdc(day, ascites),
         bili = tdc(day, bili),
         albumin = tdc(day, albumin),
         protime = tdc(day, protime), 
         alk.phos = tdc(day, alk.phos)))
