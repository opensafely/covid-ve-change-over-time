from cohortextractor import (
    StudyDefinition, 
    patients, 
    filter_codes_by_category
)

# Import codelists.py script
from codelists import *

import pandas as pd

# import the vairables for deriving JCVI groups
from grouping_variables import (
    seed,
    max_comparisons,
    jcvi_variables, 
    start_date,
    end_date,
    study_parameters,
)

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(seed)

# import recurring event functions
from recurrent_event_funs import *


study=StudyDefinition(

    default_expectations={
        "date": {"earliest": start_date, "latest": end_date},
        "rate": "uniform",
        "incidence": 0.8,
    },

    **jcvi_variables,   

    population=patients.satisfying(
        """
        has_follow_up = 1
        AND
        NOT died_before
        """,
        has_follow_up=patients.registered_with_one_practice_between(
            start_date="elig_date - 1 year",
            end_date="elig_date",
            return_expectations={"incidence": 0.90},
        ),
        died_before=patients.died_from_any_cause(
            on_or_before="elig_date - 1 day",
            returning="binary_flag",
        ),
    ),

    # as IMD and region are all defined on a given date, these can only be updated at elig_date + 8 weeks... 
    # this is the ealiest possible date for start_1 as:
    # elig_date + 6 weeks is earliest possible 2nd vax date, + 2 weeks for start of comparison 1
    
    # BMI will be derived as a sequence
    # all other covariates are "ever", so can just find the earliest ate before end_K

    #############################
    ### DEMOGRAPHIC VARIABLES ###
    #############################

    # DOB
    dob=patients.date_of_birth(
        date_format="YYYY-MM",
        return_expectations={
            "incidence": 1
        }
    ),

    # Patients in long-stay nursing and residential care
    # any time before end_K_date (return earliest date)
    longres_date=patients.with_these_clinical_events(
        longres_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_after="elig_date + 56 days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),

    imd=patients.address_as_of(
                    "elig_date + 56 days",
                    returning="index_of_multiple_deprivation",
                    round_to_nearest=100,
                    return_expectations={
                        "category": {"ratios": {c: 1/320 for c in range(100,32100,100)}}
                        }
                    ),

    # # ##########################
    # # ### CLINICAL VARIABLES ###
    # # ##########################

    # 10 most recent bmi recordings before end_K
    # bmi_1 is most recent, bmi_2 second most recent etc.
    # elig_date + 16 weeks (latest 2nd vax) + 2 weeks (latest start_1) + (k+1)*28 days (latest end_k)
    **most_recent_bmi_X(
        name="bmi",
        n=study_parameters["recur_bmi"],
        index_date=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "float": {"distribution": "normal", "mean": 28, "stddev": 8},
            "incidence": 0.80,
        },
    ),

    # chronic caridac disease
    chronic_cardiac_disease_date=patients.with_these_clinical_events(
        chronic_cardiac_disease_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
   
    # heart failure
     heart_failure_date=patients.with_these_clinical_events(
        heart_failure_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # other heart disease
    other_heart_disease_date=patients.with_these_clinical_events(
        other_heart_disease_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
   
    # diabetes
    diabetes_date=patients.with_these_clinical_events(
        diabetes_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
   
    # dialysis
    dialysis_date=patients.with_these_clinical_events(
        dialysis_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
   
    # chronic liver disease
    chronic_liver_disease_date=patients.with_these_clinical_events(
        chronic_liver_disease_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # current COPD
    current_copd_date=patients.with_these_clinical_events(
        current_copd_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # learning disability including downs synrdome and cerebal palsy
    ld_inc_ds_and_cp_date=patients.with_these_clinical_events(
        learning_disability_including_downs_syndrome_and_cerebral_palsy_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # cystic fibrosis
    cystic_fibrosis_date=patients.with_these_clinical_events(
        cystic_fibrosis_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # other respiratory conditions
    other_respiratory_date=patients.with_these_clinical_events(
        other_respiratory_conditions_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # lung cancer
    lung_cancer_date=patients.with_these_clinical_events(
        lung_cancer_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
   
    # haematological cancer
    haematological_cancer_date=patients.with_these_clinical_events(
        haematological_cancer_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
   
    # cancers excluding excluding lunch and haematological
    cancer_excl_lung_and_haem_date=patients.with_these_clinical_events(
        cancer_excluding_lung_and_haematological_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # chemo or radiotherapy
    chemo_or_radio_date=patients.with_these_clinical_events(
        chemotherapy_or_radiotherapy_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # solid organ transplant
    solid_organ_transplantation_date=patients.with_these_clinical_events(
        solid_organ_transplantation_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # bone marrow transplant
    bone_marrow_transplant_date=patients.with_these_clinical_events(
        bone_marrow_transplant_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # sickle cell disease
    sickle_cell_disease_date=patients.with_these_clinical_events(
        sickle_cell_disease_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # permanant immunosuppression
    permanant_immunosuppression_date=patients.with_these_clinical_events(
        permanent_immunosuppression_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # temporary immunosuppression
    temporary_immunosuppression_date=patients.with_these_clinical_events(
        temporary_immunosuppression_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # asplenia
    asplenia_date=patients.with_these_clinical_events(
        asplenia_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
   
    # dmards
    dmards_date=patients.with_these_clinical_events(
        dmards_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # dementia
    dementia_date=patients.with_these_clinical_events(
        dementia_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # other neuro conditions
    other_neuro_conditions_date=patients.with_these_clinical_events(
        other_neuro_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # psychosis, schizophrenia, bipolar
    psychosis_schiz_bipolar_date=patients.with_these_clinical_events(
        psychosis_schizophrenia_bipolar_affective_disease_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # end of life
    endoflife_date=patients.with_these_clinical_events(
        eol_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # midazolam
    midazolam_date=patients.with_these_medications(
        midazolam_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),

    # flu vaccine in the 5 years before elig_date
    flu_vaccine=patients.satisfying(
        """
        flu_vaccine_tpp_table>0 OR
        flu_vaccine_med>0 OR
        flu_vaccine_clinical>0
        """,
        
        flu_vaccine_tpp_table=patients.with_tpp_vaccination_record(
            target_disease_matches="INFLUENZA",
            between=["elig_date - 5 years", "elig_date"], 
            returning="binary_flag",
        ),
        
        flu_vaccine_med=patients.with_these_medications(
            flu_med_codes,
            between=["elig_date - 5 years", "elig_date"], 
            returning="binary_flag",
        ),
        flu_vaccine_clinical=patients.with_these_clinical_events(
            flu_clinical_given_codes,
            ignore_days_where_these_codes_occur=flu_clinical_not_given_codes,
            between=["elig_date - 5 years", "elig_date"], 
            returning="binary_flag",
        ),
        return_expectations={"incidence": 0.5, },
    ),

    # electronic frailty index
    # date hard-coded because there are few dates available for efi
    efi=patients.with_these_decision_support_values(
        algorithm="electronic_frailty_index",
        on_or_before="elig_date + 56 days", 
        find_last_match_in_period=True,
        returning="numeric_value",
        return_expectations={
            "float": {"distribution": "normal", "mean": 0.20, "stddev": 0.09},
            "incidence": 0.99
        },
    ),
    
    # dates of shielding codes
    shielded_0_date=patients.with_these_clinical_events(
        shield_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 56 days",
        find_last_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    **with_these_clinical_events_date_X(
        name="shielded",
        n=study_parameters["recur_shielded"],
        index_date="elig_date + 57 days",
        codelist=shield_primis,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.01,
        },
    ),
    
    # dates of non shielding codes
    nonshielded_0_date=patients.with_these_clinical_events(
        nonshield_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 56 days",
        find_last_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    **with_these_clinical_events_date_X(
        name="nonshielded",
        n=study_parameters["recur_shielded"],
        index_date="elig_date + 57 days",
        codelist=nonshield_primis,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.01,
        },
    ),

    # # # ##############
    # # # ### EVENTS ###
    # # # ##############
    
    # positive covid test
    positive_test_0_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),

    # probable covid case identified in primary care
    # Will also had 'covid in primary care', but as far as I can see it was the same as this probable definition.
    # I have used _0_ in these in case the study design changes to no longer exclude anyone with
    # previous COVID, in which case the outcome variables must become recurrent
    primary_care_covid_case_0_date=patients.with_these_clinical_events(
        covid_primary_care_probable_combined,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # suspected covid case identified in primary care
    primary_care_suspected_covid_0_date=patients.with_these_clinical_events(
        primary_care_suspected_covid_combined,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # unplanned hospital admission
    admitted_unplanned_0_date=patients.admitted_to_hospital(
        returning="date_admitted",
        on_or_before="elig_date + 56 days",
        find_first_match_in_period=True,
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        date_format="YYYY-MM-DD",
    ),
     **admitted_date_X(
        name = "admitted_unplanned",
        n = study_parameters["recur_admissions"],
        index_name = "admitted_unplanned",
        index_date = "elig_date + 57 days",
        returning="date_admitted",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
    ),
    
    discharged_unplanned_0_date=patients.admitted_to_hospital(
        returning="date_discharged",
        on_or_before="elig_date + 56 days",
        find_first_match_in_period=True,
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        date_format="YYYY-MM-DD",
    ), 
    **admitted_date_X(
        name = "discharged_unplanned",
        n = study_parameters["recur_admissions"],
        index_name = "admitted_unplanned",
        index_date = "elig_date + 57 days",
        returning="date_admitted",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
    ),
   
    # unplanned infectious hospital admission
    admitted_unplanned_infectious_0_date=patients.admitted_to_hospital(
        returning="date_admitted",
        on_or_before="elig_date + 56 days",
        find_first_match_in_period=True,
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    **admitted_date_X(
        name = "admitted_unplanned_infectious",
        n = study_parameters["recur_admissions"],
        index_name = "admitted_unplanned_infectious",
        index_date = "elig_date + 57 days",
        returning="date_admitted",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
    ),
   
    discharged_unplanned_infectious_0_date=patients.admitted_to_hospital(
        returning="date_discharged",
        on_or_before="elig_date + 56 days",
        find_first_match_in_period=True,
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        date_format="YYYY-MM-DD",
    ), 
    **admitted_date_X(
        name = "discharged_unplanned_infectious",
        n = study_parameters["recur_admissions"],
        index_name = "admitted_unplanned_infectious",
        index_date = "elig_date + 57 days",
        returning="date_discharged",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
    ),
    
    # covid hospital adamission
    covidadmitted_0_date=patients.admitted_to_hospital(
        returning="date_admitted",
        with_these_diagnoses=covid_codes,
        on_or_before=f"elig_date + {112 + 14 + (max_comparisons + 1)*28} days",
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01,
        },
    ),
    
    # covid death
    coviddeath_date=patients.with_these_codes_on_death_certificate(
        covid_codes,
        returning="date_of_death",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.02
        },
    ),
    # any death
    death_date=patients.died_from_any_cause(
        returning="date_of_death",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.02
        },
    ),
    # De-registration
    dereg_date=patients.date_deregistered_from_all_supported_practices(
        on_or_after="elig_date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date,},
            "incidence": 0.001
        }
    ),
    
)