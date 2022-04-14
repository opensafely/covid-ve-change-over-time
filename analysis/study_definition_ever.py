from cohortextractor import (
    StudyDefinition, 
    patients
)

# Import codelists.py script
from codelists import *

# import json module
import json
#study_parameters
with open("./analysis/lib/study_parameters.json") as f:
  study_parameters = json.load(f)

# define variables explicitly
start_date=study_parameters["start_date"] # start of phase 1
end_date=study_parameters["end_date"] # latest date of data

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

def covid_test_date_X(name, index_date, n, test_result):
  # covid test date (result can be "any", "positive", or "negative")
  def var_signature(name, on_or_after, test_result):
    return {
      name: patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result=test_result,
        on_or_after=on_or_after,
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        returning="date",
        date_format="YYYY-MM-DD"
      ),
    }
  variables = var_signature(
      name=f"{name}_1_date", 
      on_or_after=f"{index_date} - 42 days", 
      test_result=test_result)
  for i in range(2, n+1):
    variables.update(var_signature(
        name=f"{name}_{i}_date", 
        on_or_after=f"{name}_{i-1}_date + 1 day", 
        test_result=test_result))
  return variables


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

    ##############
    ### EVENTS ###
    ##############

    # positive tesy
    # date of earliest before start_1_date
    positive_test_pre_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="start_1_date",
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    # date of earliest after start_1_date
    positive_test_post_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_after="start_1_date + 1 day",
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),

    # probable covid case identified in primary care
    # date of earliest before start_1_date
    primary_care_covid_case_pre_date=patients.with_these_clinical_events(
        covid_primary_care_probable_combined,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="start_1_date",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    # date of earliest after start_1_date
    primary_care_covid_case_post_date=patients.with_these_clinical_events(
        covid_primary_care_probable_combined,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_after="start_1_date + 1 day",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # covid hospitalisation:
    # from ACPS in any field
    # date of earliest before start_1_date
    covidadmitted_pre_date=patients.admitted_to_hospital(
        returning="date_admitted",
        with_these_diagnoses=covid_codes,
        on_or_before="start_1_date",
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01,
        },
    ),
    # date of earliest after start_1_date
    covidadmitted_post_date=patients.admitted_to_hospital(
        returning="date_admitted",
        with_these_diagnoses=covid_codes,
        on_or_after="start_1_date + 1 day",
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01,
        },
    ),
    # positive tests for investigating hospitalisations
    # dates of up to 6 positive tests in the period starting 6 weeks before covidadmitted_post_date
    **covid_test_date_X(
        name="covidadmitted_postest",
        index_date="covidadmitted_post_date",
        n=6,
        test_result="positive"
    ),

    # from ACPS in primary field
    # date of earliest before start_1_date
    covidadmittedprimary_pre_date=patients.admitted_to_hospital(
        returning="date_admitted",
        with_these_primary_diagnoses=covid_codes,
        on_or_before="start_1_date",
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01,
        },
    ),
    # date of earliest after start_1_date
    covidadmittedprimary_post_date=patients.admitted_to_hospital(
        returning="date_admitted",
        with_these_primary_diagnoses=covid_codes,
        on_or_after="start_1_date + 1 day",
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01,
        },
    ),

    # from ECDS
    # any emergency attendance for covid
    # date of earliest before start_1_date
    covidemergency_pre_date=patients.attended_emergency_care(
        returning="date_arrived",
        with_these_diagnoses=covid_emergency,
        discharged_to=discharged_to_hospital,
        on_or_before="start_1_date",
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
    ),
    # date of earliest after start_1_date
    covidemergency_post_date=patients.attended_emergency_care(
        returning="date_arrived",
        with_these_diagnoses=covid_emergency,
        discharged_to=discharged_to_hospital,
        on_or_after="start_1_date + 1 day",
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
    ),


)