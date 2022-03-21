# Waning effectiveness of BNT162b2 and ChAdOx1 COVID-19 vaccines over six months since second dose: a cohort study using linked electronic health records

The aim of this study was to comapre rates of COVID-19 hospitalisation, COVID-19 and non-COVID-19 mortality, and infection with SARS-CoV-2, between adults fully vaccinated with the Pfizer-BioNTech BNT162b2 mRNA vaccine (BNT162b2) and the Oxford-AstraZeneca ChAdOx1 nCoV-19 AZD1222 (ChAdOx1), and those who were unvaccinated.

## Repository navigation

You can run this project via [Gitpod](https://gitpod.io) in a web browser by clicking on this badge: [![Gitpod ready-to-code](https://img.shields.io/badge/Gitpod-ready--to--code-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/opensafely/covid-ve-change-over-time)

* The [protocol is in the OpenSAFELY Google drive]()
* The pre-print is [here]()
* Analysis scripts are in the [`analysis/`](./analysis) directory (see below for details)
* Non-disclosive model outputs, including tables, figures, etc, are in `released_outputs/`
* If you are interested in how we defined our code lists, look in the [codelists folder](./codelists/)
* The [`project.yaml`](./project.yaml) defines run-order and dependencies for all the analysis scripts. **This file should *not* be edited directly**. To make changes to the yaml, edit and run the [`create-project.R`](./create-project.R) script instead.

### Analysis scripts

* [`design.R`](./analysis/design.R) creates metadata for aspects of the study design
* The following files outline how we defined our variables:
  * [`study_definition_vax.py`](analysis/study_definition_vax.py) defines the JCVI groups (see [`grouping_variables.py`](analysis/grouping_variables.py)), the vaccine vairables, and the vairables used to apply the initial eligibility criteria
  * [`study_definition_ever.py`](analysis/study_definition_ever.py) defines the covairates for the model that are either based on an "ever" diagnosis, or are not updated at the start of each comparison period
  * [`study_definition_k.py`](analysis/study_definition_k.py) is the template for definining covariates at the start of comparison *k*; see the [`create-project.R`](./create-project.R) script for how to create a script for each comparison period from this template
* [`dummy_data_vax.R`](analysis/dummy_data_vax.R) contains the script used to generate dummy data. This is used instead of the usual dummy data specified in [`study_definition_vax.py`](analysis/study_definition_vax.py), because it is then possible to impose some more useful structure in the data, such as ensuring nobody has a first dose of both the Pfizer and Astra-Zeneca vaccines. If [`study_definition_vax.py`](analysis/study_definition_vax.py) is updated, this script must also be updated to ensure variable names and types match.
* [`lib/`](./analysis/lib)
  * [`data_process_functions.R`](analysis/lib/data_process_functions.R) contains functions used in data processing
  * [`dummy_data_functions.R`](analysis/lib/dummy_data_functions.R) contains functions used to define the dummy data
  * [`data_properties.R`](analysis/lib/data_properties.R) contains functions to summarise data properties
  * [`process_covariates.R`](analysis/lib/process_covariates.R) is a script for defining covariates on `start_fu_date`
  * [`redaction_functions.R`](analysis/lib/redaction_functions.R) contains functions for redacting values <= a specified threshold
  * [`round_km.R`](analysis/lib/round_km.R) contains functions for rounding Kaplan-Meier estimates for disclosure control
* [`preprocess/`](.analysis/preprocess)

# About the OpenSAFELY framework

The OpenSAFELY framework is a Trusted Research Environment (TRE) for electronic
health records research in the NHS, with a focus on public accountability and
research quality. Read more at [OpenSAFELY.org](https://opensafely.org).

Developers and epidemiologists interested in the framework should review [the OpenSAFELY documentation](https://docs.opensafely.org)

# Licences
As standard, research projects have a MIT license. 
