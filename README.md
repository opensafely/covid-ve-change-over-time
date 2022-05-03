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

Below is a description of each of the actions in the [`project.yaml`](./project.yaml). Arguments are denoted by {arg} in the action name.

#### `dummy_data_vax`
Runs [`dummy_data_vax.R`](analysis/dummy_data_vax.R) to generate dummy data. 
This is used instead of the usual dummy data specified in [`study_definition_vax.py`](analysis/study_definition_vax.py), because it is then possible to impose some more useful structure in the data, such as ensuring nobody has a first dose of both the Pfizer and Astra-Zeneca vaccines. 
If [`study_definition_vax.py`](analysis/study_definition_vax.py) is updated, [`dummy_data_vax.R`](analysis/dummy_data_vax.R) must also be updated to ensure variable names and types match.

#### `generate_study_population`
Runs [`study_definition_vax.py`](analysis/study_definition_vax.py) to define the JCVI groups (see [`grouping_variables.py`](analysis/grouping_variables.py)), the vaccine vairables, the vairables used to apply the initial eligibility criteria, and the outcome variables (with the exception of anytest which is defined in [`study_definition_covs.py`](analysis/study_definition_covs.py)).

#### `data_input_process`
Runs [`data_input_process.R`](analysis/preprocess/data_input_process.R) to process `input_vax.feather` (generated by `generate_study_population`).

#### `data_eligible_ab`
Runs [`data_eligible_ab.R`](analysis/preprocess/data_eligible_ab.R) to apply eligibility criteria from boxes A and B in Figure 3 of the protocol.

#### `data_2nd_vax_dates`
Runs [`data_2nd_vax_dates.R`](analysis/second_vax_period/data_2nd_vax_dates.R) to identify the 28-day period in each strata during which the greatest number of second vaccine doses were received.

#### `plot_2nd_vax_dates`
Runs [`plot_2nd_vax_dates.R`](analysis/second_vax_period/plot_2nd_vax_dates.R) to plot the distribution of second vaccination period (SVP) dates in each stratum.

#### `data_eligible_cde`
Runs [`data_eligible_cde.R`](analysis/preprocess/data_eligible_cde.R) to apply eligibility criteria from boxes C, D and E in Figure 3 of the protocol.

#### `generate_covs_data`
Runs [`study_definition_covs.py`](analysis/study_definition_covs.py) to define covariates for the Cox model, which are defined at the start of the SVP.

#### `data_covariates_process`
Runs [`data_covariates_process.R`](analysis/preprocess/data_covariates_process.R) to process `input_covs.feather` (generated by `generate_covs_data`).

#### `data_min_max_fu`
Runs [`data_min_max_fu.R`](analysis/comparisons/data_min_max_fu.R) which creates a csv file containing the earliest and latest follow-up dates in each subgroup, and stores them in `data_min_max_fu.csv` for release.

#### `plot_cumulative_incidence`
Runs [`plot_cumulative_incidence.R`](analysis/subsequent_vax/plot_cumulative_incidence.R) which calculates the cumulative incidence of 3rd dose in the groups that received two doses of BNT162b2 or ChAdOx1 during the SVP, and the cumulative incidence of 1st dose in the group that remained unvaccinated at the start of comparison period 1.
The cumulative incidence is stored in `survtable_redacted.csv` for release.

#### `table1`
Runs [`plot_cumulative_incidence.R`](analysis/report/table1.R)  which derives the summary statistics for Table 1 of the manuscript.
These summary statistics are saved in `table1.csv` for release.
This action also combines the eligibility counts that have been generated from previous actions, and stores them in `eligibility_count_all.csv` for release.

#### `data_tte_process_{comparison}`
Runs [data_tte_process.R](analysis/comparisons/data_tte_process.R) which derives time to event data for each comparison.

#### `preflight_{comparison}_{subgroup}_{outcome}`
Runs [preflight.R](analysis/comparisons/preflight.R) which preprocesses the data for the Cox models, for each compairson, subgroup and outcome.

#### `apply_cox_{comparison}_{subgroup}_{outcome}`
Runs [apply_model_cox_update.R](analysis/comparisons/apply_model_cox_update.R) which fits Cox models for each compairson, subgroup and outcome.
Each script fits 12 models: unadjusted and adjusted models for each of the six comparison periods.

#### `combine_estimates`
Runs [combine_estimates.R](analysis/comparisons/combine_estimates.R) which combines estimates from all Cox models and stores them in `estimates_all.csv` for release.
This script also combines tables of event counts from [data_tte_process.R](analysis/comparisons/data_tte_process.R) and stores them in `event_counts_all.csv` for release.

#### `plot_check`
Runs [plot_cox_check.R](analysis/comparisons/plot_cox_check.R) which plots all estimates such that they can be checked prior to release of `estimates_all.csv`.

### `analysis/functions`
The following scripts contain functions which are used throughout this project:
* [`data_process_functions.R`](analysis/lib/data_process_functions.R) contains functions used in data processing
* [`dummy_data_functions.R`](analysis/lib/dummy_data_functions.R) contains functions used to define the dummy data
* [`data_properties.R`](analysis/lib/data_properties.R) contains functions to summarise data properties
* [`redaction_functions.R`](analysis/lib/redaction_functions.R) contains functions for redacting values <= a specified threshold
* [`round_km.R`](analysis/lib/round_km.R) contains functions for rounding Kaplan-Meier estimates for disclosure control


### `analysis/post_release`
This folder contains scripts that apply further analyses or post-processing of files that have been released from OpenSAFELY (and therefore are not run via the [`project.yaml`](./project.yaml).
The run-order of these scripts is given in [instructions.R](analysis/post_release/instructions.R).

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
