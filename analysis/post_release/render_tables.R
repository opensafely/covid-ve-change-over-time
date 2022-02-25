# render tables

rmarkdown::render(
  here::here("analysis","post_release", "table1_process.Rmd"),
  output_file = here::here("release20220221", "table1_process.docx"))

rmarkdown::render(
  here::here("analysis","post_release", "table_covariate_estimates_bind.Rmd"),
  output_file = here::here("release20220221", "table_covariate_estimates_bind.docx"))

