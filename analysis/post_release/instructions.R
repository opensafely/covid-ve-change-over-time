# define release folder
release_folder <- here::here("output", "release_objects")


### Flow chart
source(here::here("analysis", "post_release", "flow.R")) 

### Table 1
source(here::here("analysis", "post_release", "table1_process"))

rmarkdown::render(
  here::here("analysis","post_release", "table1_process.Rmd"),
  knit_root_dir = release_folder,
    output_file = here::here(release_folder, "table1_process.docx"))

### subsequent vax
source(here::here("analysis", "subsequent_vax", "plot_cumulative_incidence.R"))

### Metaregression

# Run waning_metareg.do in Stata

source(here::here("analysis", "post_release", "data_metareg_process.R"))

source(here::here("analysis", "post_release", "data_metareg_k.R"))

### Results plots

source(here::here("analysis", "post_release", "plot_cox_all.R"))

### Appendix