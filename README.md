# Waning effectiveness of BNT162b2 and ChAdOx1 COVID-19 vaccines over six months since second dose: a cohort study using linked electronic health records

The aim of this study was to comapre rates of COVID-19 hospitalisation, COVID-19 and non-COVID-19 mortality, and infection with SARS-CoV-2, between adults fully vaccinated with the Pfizer-BioNTech BNT162b2 mRNA vaccine (BNT162b2) and the Oxford-AstraZeneca ChAdOx1 nCoV-19 AZD1222 (ChAdOx1), and those who were unvaccinated.

## Repository navigation

You can run this project via [Gitpod](https://gitpod.io) in a web browser by clicking on this badge: [![Gitpod ready-to-code](https://img.shields.io/badge/Gitpod-ready--to--code-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/opensafely/covid-ve-change-over-time)

* The [protocol is in the OpenSAFELY Google drive]()
* The pre-print is [here](https://www.medrxiv.org/content/10.1101/2022.03.23.22272804v1)
* Analysis scripts are in the [`analysis/`](./analysis) directory (see below for details)
* Non-disclosive model outputs, including tables, figures, etc, are in `released_outputs/`
* If you are interested in how we defined our code lists, look in the [codelists folder](./codelists/)

### project.yaml
The [`project.yaml`](./project.yaml) defines run-order and dependencies for all the analysis scripts. 
**This file should *not* be edited directly**. To make changes to the yaml, edit and run the [`create-project.R`](./create-project.R) script instead.
The [`create-project.R`](./create-project.R) also creates metadata files and stores them in the [`analysis/lib`](./analysis/lib) folder.
These metadata files are annotated in [`create-project.R`](./create-project.R).
There is no need to run [`create_project.yaml`](./create_project.R) if you are simply cloning this repo.

### `analysis/`

* The following files outline how we defined our variables:
  * [`study_definition_vax.py`](analysis/study_definition_vax.py) defines the JCVI groups (see [`grouping_variables.py`](analysis/grouping_variables.py)), the vaccine vairables, and the vairables used to apply the initial eligibility criteria
  * [`study_definition_ever.py`](analysis/study_definition_ever.py) defines the covairates for the model that are either based on an "ever" diagnosis, or are not updated at the start of each comparison period
  * [`study_definition_k.py`](analysis/study_definition_k.py) is the template for definining covariates at the start of comparison *k*; see the [`create-project.R`](./create-project.R) script for how to create a script for each comparison period from this template
* [`dummy_data_vax.R`](analysis/dummy_data_vax.R) contains the script used to generate dummy data. This is used instead of the usual dummy data specified in [`study_definition_vax.py`](analysis/study_definition_vax.py), because it is then possible to impose some more useful structure in the data, such as ensuring nobody has a first dose of both the Pfizer and Astra-Zeneca vaccines. If [`study_definition_vax.py`](analysis/study_definition_vax.py) is updated, this script must also be updated to ensure variable names and types match.
* [`lib/`](./analysis/lib)
  * [`data_process_functions.R`](analysis/lib/data_process_functions.R) contains functions used in data processing
  * [`dummy_data_functions.R`](analysis/lib/dummy_data_functions.R) contains functions used to define the dummy data
  * [`data_properties.R`](analysis/lib/data_properties.R) contains functions to summarise data properties
  * [`redaction_functions.R`](analysis/lib/redaction_functions.R) contains functions for redacting values <= a specified threshold
  * [`round_km.R`](analysis/lib/round_km.R) contains functions for rounding Kaplan-Meier estimates for disclosure control
* The scripts in [`preprocess/`](.analysis/preprocess), [`second_vax_period/`](./analysis/second_vax_period), [`subsequent_vax/`](./analysis/subsequent_vax), [`comparisons/`](./analysis/comparisons), and [`report/`](./analysis/report) carry out preprocessing, analysis and postprocessing, and are applied in the [`project.yaml`](./project.yaml)
* [`post_release/`](./analysis/post_release) contains scripts that apply postprocessing to the released results:
  * [`flow.R`](analysis/post_release/flow.R) creates a table of the flow of individuals into the study
  * [`table1_process.Rmd`](analysis/post_release/table1_process.Rmd) creates Table 1 for the manuscript, which can be rendered by running [`render_tables.R`](analysis/post_release/render_tables.R)
  * [`data_metareg_process.R`](analysis/post_release/data_metareg_process.R) prepares the data for metaregression
  * [`waning_metareg.do`](analysis/post_release/waning_metareg.do) applies the metaregression
  * [`data_metareg_k.R`](analysis/post_release/data_metareg_k.R) applies some postprocessing to the metaregression results
  * [`plot_cox_all.R`](analysis/post_release/plot_cox_all.R) plots the results of the Cox regression and metaregression

### `manuscript/`

* [`manuscript_text.Rmd`](manuscript/manuscript_text.Rmd) contains some paragraphs from the manuscript that have a lot of references to the effect estimates **note that this script needs updating with updated wording**
* [`appendix.Rmd`](manuscript/appendix.Rmd) contains the supplementary material for the paper; to render to a PDF requires [`preamble.tex`](manuscript/preamble.tex) and [`preface.tex`](manuscript/preface.tex)

# About the OpenSAFELY framework

The OpenSAFELY framework is a Trusted Research Environment (TRE) for electronic
health records research in the NHS, with a focus on public accountability and
research quality. Read more at [OpenSAFELY.org](https://opensafely.org).

Developers and epidemiologists interested in the framework should review [the OpenSAFELY documentation](https://docs.opensafely.org)

# Licences
As standard, research projects have a MIT license. 
