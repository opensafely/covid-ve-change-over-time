# render tables

rmarkdown::render(
  here::here("analysis","post_release", "table1_process.Rmd"),
  output_file = here::here("release20220221", "table1_process.docx"))

rmarkdown::render(
  here::here("manuscript", "table_covariate_estimates.Rmd"),
  output_file = here::here("release20220221", "table_covariate_estimates.docx"))

