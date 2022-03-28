from cohortextractor import (
    StudyDefinition, 
    patients
)

# Import codelists.py script
from codelists import *

# import json module
import json
#study_parameters
with open("./output/lib/study_parameters.json") as f:
  study_parameters = json.load(f)

# define variables explicitly
start_date=study_parameters["start_date"] # start of phase 1
end_date=study_parameters["end_date"] # latest date of data

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

###
study=StudyDefinition(

    default_expectations={
        "date": {"earliest": start_date, "latest": end_date},
        "rate": "uniform",
        "incidence": 0.8,
    },  

    population=patients.all(),

    # min elig date within subgroup
    min_elig_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning='min_elig_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),
    # start date of first comparison
    start_1_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning='start_1_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),

    # age on start date
    age=patients.age_as_of(
        "start_1_date",
        return_expectations={
            "rate" : "universal",
            "int" : {"distribution" : "population_ages"}
            }
        ),

    ### covid tests as covariates
    # during unvaccinated time (from when tests widely availabe to elig_date)
    test_hist_n=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="any",
        between=["2020-05-18", "min_elig_date - 1 day"], # day before 1st vaccine eligibility date
        restrict_to_earliest_specimen_date=False,
        returning="number_of_matches_in_period",
        return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	    ),

    #### ever covariates
    # Asthma Diagnosis code ever
    astdx_date=patients.with_these_clinical_events(
        ast_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
    ),

    # Chronic Respiratory Disease other than asthma
    resp_date=patients.with_these_clinical_events(
        resp_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
    ),

    # Chronic Neurological Disease including Significant Learning Disorder
    cns_date=patients.with_these_clinical_events(
        cns_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
    ),

    # Diabetes
    diab_date=patients.with_these_clinical_events(
        diab_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
        ),

    # Severe mental illness
    sev_mental_date=patients.with_these_clinical_events(
        sev_mental_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
        ),

    # Chronic heart disease codes
    chd_date=patients.with_these_clinical_events(
        chd_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
    ),

    # Chronic Liver disease codes
    cld_date=patients.with_these_clinical_events(
        cld_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
    ),

    # Immunosuppression diagnosis codes
    immdx_date=patients.with_these_clinical_events(
        immdx_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
        ),

    # Asplenia or Dysfunction of the Spleen codes
    spln_date=patients.with_these_clinical_events(
        spln_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
    ),

    # Learning Disability
    learndis_date=patients.with_these_clinical_events(
        learndis_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
    ),

    # Patients in long-stay nursing and residential care before start period k
    longres_date=patients.with_these_clinical_events(
        longres_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=end_date,
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.02},
    ),

    # flu vaccine in the 5 years before start_k_date
    flu_vaccine=patients.satisfying(
        """
        flu_vaccine_tpp_table>0 OR
        flu_vaccine_med>0 OR
        flu_vaccine_clinical>0
        """,
        
        flu_vaccine_tpp_table=patients.with_tpp_vaccination_record(
            target_disease_matches="INFLUENZA",
            between=["start_1_date - 5 years", "start_1_date"], 
            returning="binary_flag",
        ),
        
        flu_vaccine_med=patients.with_these_medications(
            flu_med_codes,
            between=["start_1_date - 5 years", "start_1_date"], 
            returning="binary_flag",
        ),
        flu_vaccine_clinical=patients.with_these_clinical_events(
            flu_clinical_given_codes,
            ignore_days_where_these_codes_occur=flu_clinical_not_given_codes,
            between=["start_1_date - 5 years", "start_1_date"], 
            returning="binary_flag",
        ),
        return_expectations={"incidence": 0.5, },
    ),

  # any emergency attendance for covid
  covidemergency_0_date=patients.attended_emergency_care(
    returning="date_arrived",
    on_or_before=end_date,
    with_these_diagnoses=covid_emergency,
    discharged_to=discharged_to_hospital,
    date_format="YYYY-MM-DD",
    find_first_match_in_period=True,
  ),
  

)