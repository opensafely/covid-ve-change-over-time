# Investigating the distribution of positive tests around COVID-19 hospitalisations

# create output folder
fs::dir_create(here::here("output", "eda"))

# render report 
rmarkdown::render(
  "analysis/eda/covidadmitted_postest.Rmd",
  output_file="covidadmitted_postest",
  output_dir="output/eda")
