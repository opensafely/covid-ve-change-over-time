from cohortextractor import (
    StudyDefinition, 
    patients, 
    filter_codes_by_category
)

# Import codelists.py script
from codelists import *

# import the vairables for deriving JCVI groups
from grouping_variables import (
    jcvi_variables, 
    start_date,
    end_date,
)

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

    #############################
    ### DEMOGRAPHIC VARIABLES ###
    #############################

    # dob=patients.date_of_birth(
    #     date_format="YYYY-MM",
    #     return_expectations={
    #         "incidence": 1
    #     }
    # ),
    
     # hscworker=patients.with_healthcare_worker_flag_on_covid_vaccine_record(
     #     returning="binary_flag",
     #     return_expectations={"incidence": 0.01},
     #     ),

    # ETHNICITY IN 6 CATEGORIES
    # ethnicity
    ethnicity_6=patients.with_these_clinical_events(
        eth2001_primis,
        returning="category",
        find_last_match_in_period=True,
        on_or_before="elig_date + 42 days",
        return_expectations={
            "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
            "incidence": 0.75,
        },
    ),

    ethnicity_6_sus = patients.with_ethnicity_from_sus(
        returning = "group_6",  
        use_most_frequent_code = True,
        return_expectations = {
            "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
            "incidence": 0.8,
        },
    ),

    # IMD - quintile
    imd=patients.categorised_as(
        {
            "0": "DEFAULT",
            "1": """index_of_multiple_deprivation >=1 AND index_of_multiple_deprivation < 32844*1/5""",
            "2": """index_of_multiple_deprivation >= 32844*1/5 AND index_of_multiple_deprivation < 32844*2/5""",
            "3": """index_of_multiple_deprivation >= 32844*2/5 AND index_of_multiple_deprivation < 32844*3/5""",
            "4": """index_of_multiple_deprivation >= 32844*3/5 AND index_of_multiple_deprivation < 32844*4/5""",
            "5": """index_of_multiple_deprivation >= 32844*4/5 """,
        },
        index_of_multiple_deprivation=patients.address_as_of(
            "elig_date + 42 days",
            returning="index_of_multiple_deprivation",
            round_to_nearest=100,
        ),
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {
                    "0": 0.01,
                    "1": 0.20,
                    "2": 0.20,
                    "3": 0.20,
                    "4": 0.20,
                    "5": 0.19,
                }
            },
        },
    ),

    # region - NHS England 9 regions
    region=patients.registered_practice_as_of(
        "elig_date + 42 days",
        returning="nuts1_region_name",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {
                    "North East": 0.1,
                    "North West": 0.1,
                    "Yorkshire and The Humber": 0.1,
                    "East Midlands": 0.1,
                    "West Midlands": 0.1,
                    "East": 0.1,
                    "London": 0.2,
                    "South West": 0.1,
                    "South East": 0.1
                },
            },
            "incidence": 0.99
        },
    ),

    # Patients in long-stay nursing and residential care
    longres_0_date=patients.with_these_clinical_events(
        longres_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 42 days",
        find_last_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    # **with_these_clinical_events_date_X(
    #     name="longres",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=longres_primis,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),

    ##########################
    ### CLINICAL VARIABLES ###
    ##########################

    # BMI
    bmi_0=patients.most_recent_bmi(
        on_or_before="elig_date + 42 days",
        minimum_age_at_measurement=16,
        # on_most_recent_day_of_measurement=True, # returning an error for some reason
        include_measurement_date=True,
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "float": {"distribution": "normal", "mean": 28, "stddev": 8},
            "incidence": 0.80,
        },
    ),
    # **most_recent_bmi_X(
    #     name="bmi",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "float": {"distribution": "normal", "mean": 28, "stddev": 8},
    #         "incidence": 0.80,
    #     },
    # ),
    # 
    # # chronic caridac disease
    # chronic_cardiac_disease_0_date=patients.with_these_clinical_events(
    #     chronic_cardiac_disease_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="chronic_cardiac_disease",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=chronic_cardiac_disease_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # heart failure
    #  heart_failure_0_date=patients.with_these_clinical_events(
    #     heart_failure_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="heart_failure",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=heart_failure_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # other heart disease
    # other_heart_disease_0_date=patients.with_these_clinical_events(
    #     other_heart_disease_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="other_heart_disease",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=other_heart_disease_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # diabetes
    # diabetes_0_date=patients.with_these_clinical_events(
    #     diabetes_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="diabetes",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=diabetes_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # dialysis
    # dialysis_0_date=patients.with_these_clinical_events(
    #     dialysis_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="dialysis",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=dialysis_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # chronic liver disease
    # chronic_liver_disease_0_date=patients.with_these_clinical_events(
    #     chronic_liver_disease_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="chronic_liver_disease",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=chronic_liver_disease_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # current COPD
    # current_copd_0_date=patients.with_these_clinical_events(
    #     current_copd_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="current_copd",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=current_copd_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # learning disability including downs synrdome and cerebal palsy
    # ld_inc_ds_and_cp_0_date=patients.with_these_clinical_events(
    #     learning_disability_including_downs_syndrome_and_cerebral_palsy_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="ld_inc_ds_and_cp",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=learning_disability_including_downs_syndrome_and_cerebral_palsy_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # cystic fibrosis
    # cystic_fibrosis_0_date=patients.with_these_clinical_events(
    #     cystic_fibrosis_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="cystic_fibrosis",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=cystic_fibrosis_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # other respiratory conditions
    # other_respiratory_0_date=patients.with_these_clinical_events(
    #     other_respiratory_conditions_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="other_respiratory",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=other_respiratory_conditions_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # lung cancer
    # lung_cancer_0_date=patients.with_these_clinical_events(
    #     lung_cancer_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="lung_cancer",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=lung_cancer_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # haematological cancer
    # haematological_cancer_0_date=patients.with_these_clinical_events(
    #     haematological_cancer_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="haematological_cancer",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=haematological_cancer_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # cancers excluding excluding lunch and haematological
    # cancer_excl_lung_and_haem_0_date=patients.with_these_clinical_events(
    #     cancer_excluding_lung_and_haematological_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="cancer_excl_lung_and_haem",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=cancer_excluding_lung_and_haematological_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # chemo or radiotherapy
    # chemo_or_radio_0_date=patients.with_these_clinical_events(
    #     chemotherapy_or_radiotherapy_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="chemo_or_radio",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=chemotherapy_or_radiotherapy_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # solid organ transplant
    # solid_organ_transplantation_0_date=patients.with_these_clinical_events(
    #     solid_organ_transplantation_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="solid_organ_transplantation",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=solid_organ_transplantation_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # bone marrow transplant
    # bone_marrow_transplant_0_date=patients.with_these_clinical_events(
    #     bone_marrow_transplant_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="bone_marrow_transplant",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=bone_marrow_transplant_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # sickle cell disease
    # sickle_cell_disease_0_date=patients.with_these_clinical_events(
    #     sickle_cell_disease_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="sickle_cell_disease",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=sickle_cell_disease_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # permanant immunosuppression
    # permanant_immunosuppression_0_date=patients.with_these_clinical_events(
    #     permanent_immunosuppression_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="permanant_immunosuppression",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=permanent_immunosuppression_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # temporary immunosuppression
    # temporary_immunosuppression_0_date=patients.with_these_clinical_events(
    #     temporary_immunosuppression_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="temporary_immunosuppression",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=temporary_immunosuppression_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # asplenia
    # asplenia_0_date=patients.with_these_clinical_events(
    #     asplenia_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="asplenia",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=asplenia_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # dmards
    # dmards_0_date=patients.with_these_clinical_events(
    #     dmards_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="dmards",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=dmards_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # dementia
    # dementia_0_date=patients.with_these_clinical_events(
    #     dementia_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="dementia",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=dementia_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # other neuro conditions
    # other_neuro_conditions_0_date=patients.with_these_clinical_events(
    #     other_neuro_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="other_neuro_conditions",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=other_neuro_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # psychosis, schizophrenia, bipolar
    # psychosis_schiz_bipolar_0_date=patients.with_these_clinical_events(
    #     psychosis_schizophrenia_bipolar_affective_disease_codes,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="psychosis_schiz_bipolar",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=psychosis_schizophrenia_bipolar_affective_disease_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),

    # end of life
    endoflife_0_date=patients.with_these_clinical_events(
        eol_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 42 days",
        find_last_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    # **with_these_clinical_events_date_X(
    #     name="endoflife",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=eol_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),

    # midazolam
    midazolam_0_date=patients.with_these_medications(
        midazolam_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 42 days",
        find_last_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    # **with_these_medications_date_X(
    #     name="midazolam",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=midazolam_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),

    # # flu vaccine in the last 5 years
    # flu_vaccine=patients.satisfying(
    #     """
    #     flu_vaccine_tpp_table>0 OR
    #     flu_vaccine_med>0 OR
    #     flu_vaccine_clinical>0
    #     """,
    #     
    #     flu_vaccine_tpp_table=patients.with_tpp_vaccination_record(
    #         target_disease_matches="INFLUENZA",
    #         between=["elig_date - 5 years", "elig_date"], 
    #         returning="binary_flag",
    #     ),
    #     
    #     flu_vaccine_med=patients.with_these_medications(
    #         flu_med_codes,
    #         between=["elig_date - 5 years", "elig_date"], 
    #         returning="binary_flag",
    #     ),
    #     flu_vaccine_clinical=patients.with_these_clinical_events(
    #         flu_clinical_given_codes,
    #         ignore_days_where_these_codes_occur=flu_clinical_not_given_codes,
    #         between=["elig_date - 5 years", "elig_date"], 
    #         returning="binary_flag",
    #     ),
    #     return_expectations={"incidence": 0.5, },
    # ),
    # 
    # # electronic frailty index
    # # date currently hard-coded because there are no other dates available for efi
    # # check that this still the case, otherwise change to recurrent variable like those above
    # efi=patients.with_these_decision_support_values(
    #     algorithm = "electronic_frailty_index",
    #     on_or_before = "2020-12-08", 
    #     find_last_match_in_period = True,
    #     returning="numeric_value",
    #     return_expectations={
    #         #"category": {"ratios": {0.1: 0.25, 0.15: 0.25, 0.30: 0.25, 0.5: 0.25}},
    #         "float": {"distribution": "normal", "mean": 0.20, "stddev": 0.09},
    #         "incidence": 0.99
    #     },
    # ),
    # 
    # # dates of shielding codes
    # shielded_0_date=patients.with_these_clinical_events(
    #     shield_primis,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="shielded",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=shield_primis,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),
    # 
    # # dates of non shielding codes
    # nonshielded_0_date=patients.with_these_clinical_events(
    #     nonshield_primis,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **with_these_clinical_events_date_X(
    #     name="nonshielded",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=nonshield_primis,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.01,
    #     },
    # ),

    ######################
    ### COVID VACCINES ###
    ######################

    ## any covid vaccination, identified by target disease
    covid_vax_disease_1_date=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        on_or_after=start_date,
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),
    covid_vax_disease_2_date=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        on_or_after="covid_vax_disease_1_date + 1 day",
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),
    covid_vax_disease_3_date=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        on_or_after="covid_vax_disease_2_date + 1 day",
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
         return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),

    # Pfizer BioNTech - first record of a pfizer vaccine 
    # NB *** may be patient's first COVID vaccine dose or their second if mixed types are given ***
       
    covid_vax_pfizer_1_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
        on_or_after=start_date,  
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
         return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ), 
    covid_vax_pfizer_2_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
        on_or_after="covid_vax_pfizer_1_date + 1 day",  
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),
    covid_vax_pfizer_3_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
         on_or_after="covid_vax_pfizer_2_date + 1 day",  
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
         return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),
    
    ## Oxford AZ - first record of an Oxford AZ vaccine 
    # NB *** may be patient's first COVID vaccine dose or their second if mixed types are given ***
    covid_vax_az_1_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV",
        on_or_after=start_date,
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),
    covid_vax_az_2_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV",
        on_or_after="covid_vax_az_1_date + 1 day",  
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),
    covid_vax_az_3_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV",
        on_or_after="covid_vax_az_2_date + 1 day",  
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
         return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),
    
    ## Moderna - first record of moderna vaccine
    ## NB *** may be patient's first COVID vaccine dose or their second if mixed types are given ***
    covid_vax_moderna_1_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
        on_or_after=start_date,
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
         return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),            
    covid_vax_moderna_2_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
        on_or_after="covid_vax_moderna_1_date + 1 day",  
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
         return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),
    covid_vax_moderna_3_date=patients.with_tpp_vaccination_record(
        product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
        on_or_after="covid_vax_moderna_2_date + 1 day",  
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {
                "earliest": start_date,  
                "latest": end_date,
            },
            "incidence": 0.5
        },
    ),

    ##############
    ### EVENTS ###
    ##############

    # covid test
    # **covid_test_date_X(
    #     name="covid_test",
    #     index_date="elig_date + 43 days",
    #     n=6,
    #     test_result="any",
    #     return_expectations = {
    #         "date": {"earliest": start_date,  "latest" : end_date},
    #         "rate": "exponential_increase",
    #     },
    # ),
    
    # positive covid test
    positive_test_0_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 42 days",
        find_last_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    # **covid_test_date_X(
    #     name="positive_test",
    #     index_date="elig_date + 43 days",
    #     n=6,
    #     test_result="positive",
    #     return_expectations={
    #         "date": {"earliest": start_date,  "latest" : end_date},
    #         "rate": "exponential_increase",
    #     },
    # ),

    # probable covid case identified in primary care
    # Will also had 'covid in primary care', but as far as I can see it was the same as this probable definition.
    primary_care_covid_case_0_date=patients.with_these_clinical_events(
        covid_primary_care_probable_combined,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 42 days",
        find_last_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # **with_these_clinical_events_date_X(
    #     name="primary_care_covid_case",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=covid_primary_care_probable_combined,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    # 
    # suspected covid case identified in primary care
    primary_care_suspected_covid_0_date=patients.with_these_clinical_events(
        primary_care_suspected_covid_combined,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 42 days",
        find_last_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    # **with_these_clinical_events_date_X(
    #     name="primary_care_suspected_covid",
    #     n=6,
    #     index_date="elig_date + 43 days",
    #     codelist=primary_care_suspected_covid_combined,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    
    # emergency attendance
    # emergency_attendance_0_date=patients.attended_emergency_care(
    #     returning="date_arrived",
    #     on_or_before="elig_date + 42 days",
    #     find_last_match_in_period=True,
    #     date_format="YYYY-MM-DD",
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "exponential_increase",
    #         "incidence": 0.01
    #     },
    # ),
    # **emergency_attendance_date_X(
    #     name = "emergency",
    #     n = 6,
    #     index_date = "elig_date + 43 days",
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    
    # unplanned hospital admission
    # admitted_unplanned_0_date=patients.admitted_to_hospital(
    #     returning="date_admitted",
    #     on_or_before="elig_date + 42 days",
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     with_patient_classification = ["1"],
    #     date_format="YYYY-MM-DD",
    #     find_first_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    # **admitted_date_X(
    #     name = "admitted_unplanned",
    #     n = 6,
    #     index_name = "admitted_unplanned",
    #     index_date = "elig_date + 43 days",
    #     returning="date_admitted",
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     with_patient_classification = ["1"],
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    # discharged_unplanned_0_date=patients.admitted_to_hospital(
    #     returning="date_discharged",
    #     on_or_before="elig_date + 42 days",
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     with_patient_classification = ["1"],
    #     date_format="YYYY-MM-DD",
    #     find_first_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ), 
    # **admitted_date_X(
    #     name = "discharged_unplanned",
    #     n = 6,
    #     index_name = "admitted_unplanned",
    #     index_date = "elig_date + 43 days",
    #     returning="date_discharged",
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     with_patient_classification = ["1"],
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    # 
    # # unplanned infectious hospital admission
    # admitted_unplanned_infectious_0_date=patients.admitted_to_hospital(
    #     returning="date_admitted",
    #     on_or_before="elig_date + 42 days",
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     with_patient_classification = ["1"],
    #     with_these_diagnoses = ICD10_I_codes,
    #     date_format="YYYY-MM-DD",
    #     find_first_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    # **admitted_date_X(
    #     name = "admitted_unplanned_infectious",
    #     n = 6,
    #     index_name = "admitted_unplanned_infectious",
    #     index_date = "elig_date + 43 days",
    #     returning="date_admitted",
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     with_patient_classification = ["1"],
    #     with_these_diagnoses = ICD10_I_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    # discharged_unplanned_infectious_0_date=patients.admitted_to_hospital(
    #     returning="date_discharged",
    #     on_or_before="elig_date + 42 days",
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     with_patient_classification = ["1"],
    #     with_these_diagnoses = ICD10_I_codes,
    #     date_format="YYYY-MM-DD",
    #     find_first_match_in_period=True,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ), 
    # **admitted_date_X(
    #     name = "discharged_unplanned_infectious",
    #     n = 6,
    #     index_name = "admitted_unplanned_infectious",
    #     index_date = "elig_date + 43 days",
    #     returning="date_discharged",
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     with_patient_classification = ["1"],
    #     with_these_diagnoses = ICD10_I_codes,
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    
    # covid hospital adamission
    covidadmitted_0_date=patients.admitted_to_hospital(
        returning="date_admitted",
        with_these_diagnoses=covid_codes,
        on_or_before="elig_date + 42 days",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01,
        },
    ),
    # **admitted_date_X(
    #     name = "covidadmitted",
    #     n = 6,
    #     index_name = "covidadmitted",
    #     index_date = "elig_date + 42 days",
    #     returning="date_admitted",
    #     with_these_diagnoses=covid_codes,
    #     with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.05,
    #     },
    # ),
    
    # # covid death
    # coviddeath_date=patients.with_these_codes_on_death_certificate(
    #     covid_codes,
    #     returning="date_of_death",
    #     date_format="YYYY-MM-DD",
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.02
    #     },
    # ),
    # # any death
    # death_date=patients.died_from_any_cause(
    #     returning="date_of_death",
    #     date_format="YYYY-MM-DD",
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date},
    #         "rate": "uniform",
    #         "incidence": 0.02
    #     },
    # ),
    # 
    # # De-registration
    # dereg_date=patients.date_deregistered_from_all_supported_practices(
    #     on_or_after="elig_date",
    #     date_format="YYYY-MM-DD",
    #     return_expectations={
    #         "date": {"earliest": start_date, "latest": end_date,},
    #         "incidence": 0.001
    #     }
    # ),

)
