tribble(
  ~type, ~name, ~description, ~encoding, ~date_defined,
  #
  "Exclusion criteria",
  "Evidence of prior COVID-19 infection",
  "Date of code corresponding to COVID-19 hospitalisation, positive SARS-CoV-2 test, or probable COVID-19 code in GP record.",
  "Date",
  "First occurence",
  #
  "Exclusion criteria",
  "End-of-life care pathway initiated",
  "Date of code corresponding to end-of-life care or midazolam injection used in end-of-life care in patient’s record.",
  "Date",
  "First occurence prior to SVP",
  #
  "Exclusion criteria",
  "Resident in care home",
  "Date of code corresponding to care home in patient’s record.",
  "Date",
  "First occurence prior to SVP",
  #
  "Strata",
  "JCVI group",
  "Priority groups defined by the JCVI expert advisory group",
  "2-12",
  "Age defined on 31 March 2021 for individuals eligible in phase 1 and 1 July 2021 for individuals eligible in phase 2.
Individuals defined as “clinically extremely vulnerable” on 18 Januray 2021 and “at risk” on 15 February 2021.
",
  #
  "Strata",
  "Eligibility date",
  "Date on which JCVI group (or age range within group) became eligible for 1st dose of COVID-19 vaccination.",
  "Date",
  "-",
  #
  "Strata",
  "Region",
  "NHS region based on practice address.",
  "East of England; Midlands; London; North East; Yorkshire; North West; South East; South West.",
  "Eligibility date + 6 weeks.",
  #
  "Demographic",
  "Sex",
  "-",
  "M; F.",
  "-",
  #
  "Demographic",
  "Age",
  "Age in whole years.",
  "Integer",
  "start_{1i}",
  #
  "Demographic",
  "Ethnicity",
  "From primary care records or SUS if missing from primary care.",
  "Black; Mixed; South Asian; White; Other.",
  "Eligibility date + 6 weeks ",
  #
  "Demographic",
  "Index of Multiple Deprivation (IMD)",
  "IMD quintile based on patient address",
  "1; 2; 3; 4; 5.",
  "Eligibility date + 6 weeks ",
  #
  "Clinical",
  "Body Mass Index (BMI)",
  "Based on numeric BMI data from primary care records.",
  "Under 30kg/m\\textsuperscipt{2} or not recorded; 30-34.9; 35-39.9; 40+ kg/m\\textsuperscipt{2}.",
  "start_{ki}",
  #
  "Clinical",
  "Chronic heart disease; Chronic respiratory disease; Chronic liver disease; Chronic kidney disease; Diabetes; Chronic neurological disease; Immunosuppressed; Learning disability including Down’s syndrome; Serious mental illness.",
  "Any code in primary care record before date defined.",
  "0; 1.",
  "start_{ki}",
  #
  "Clinical",
  "Multimorbidity",
  "Number of comorbid conditions in different organ systems.",
  "0; 1; 2+",
  "start_{ki}",
  #
  "Clinical",
  "Influenza vaccination",
  "Any record in past 5 years.",
  "0; 1.",
  "Eligibility date",
  #
  "Clinical",
  "Clinically extremely vulnerable",
  "UK government COVID-19 shielding criteria met.",
  "0; 1.",
  "start_{ki}",
  #
  "Clinical",
  "Pregnancy",
  "Pregnancy recorded in 36 weeks prior to date and no delivery code recorded more recently that pregnancy code.",
  "0; 1.",
  "start_{ki}",
  #
  "Clinical",
  "Care home",
  "Code corresponding to care home in patient’s record before evaluation date.",
  "0; 1.",
  "start_{ki}",
  #
  "Clinical",
  "Housebound",
  "Code indicating that the individual is housebound.",
  "0; 1.",
  "start_{ki}",
  #
  "Clinical",
  "Number of COVID-19 tests.",
  "SARS-CoV-2 tests were identified using SGSS records and based on swab date. Both polymerase chain reaction (PCR) and lateral flow tests will be included, without differentiation between symptomatic and asymptomatic infection.",
  "Integer",
  "Between 18 May 2020 (when widespread testing became available in England) and the earliest date of eligibility for first dose in the subgroup."
  )
