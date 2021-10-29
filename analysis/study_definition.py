from cohortextractor import (
    StudyDefinition, 
    patients, 
    filter_codes_by_category
)

# Import codelists.py script
from codelists import *

# import the vairables for deriving JCVI groups
from variables import (
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
        AND 
        NOT hscworker
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
        hscworker=patients.with_healthcare_worker_flag_on_covid_vaccine_record(
            returning="binary_flag",
            return_expectations={"incidence": 0.01},
        ),
    ),

    ####################
    ### DEMOGRAPHICS ###
    ####################

    # ETHNICITY IN 6 CATEGORIES
    # ethnicity
    ethnicity_6=patients.with_these_clinical_events(
        eth2001,
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
    imd=patients.address_as_of(
        "elig_date + 42 days",
        returning="index_of_multiple_deprivation",
        round_to_nearest=100,
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {c: 1/320 for c in range(100, 32100, 100)}},
            "incidence": 1,
            }    
    ),

    # region - NHS England 9 regions
    region=patients.registered_practice_as_of(
        "elig_date + 42 days",
        returning = "nuts1_region_name",
        return_expectations = {
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
        },
    ),

    # Patients in long-stay nursing and residential care
    # to capture patients entering/leaving long-stay nursing and residential care, maybe not necessary
    # most recent date before elig_date + 6 weeks
    longres_date_before=patients.with_these_clinical_events(
        longres,
        returning="date",
        on_or_before="elig_date + 42 days",
        find_last_match_in_period=True,
        return_expectations={"incidence": 0.01},
    ),
    # earliest date after elig_date + 6 weeks
    longres_date_after=patients.with_these_clinical_events(
        longres,
        returning="date",
        on_or_after="elig_date + 42 days",
        find_first_match_in_period=True,
        return_expectations={"incidence": 0.01},
    ),

    ######################
    ### COVID VACCINES ###
    ######################

    ## any covid vaccination, identified by target disease
    covid_vax_disease_1_date = patients.with_tpp_vaccination_record(
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
        product_name_matches="COVID-19 mRNA (nucleoside modified) Vaccine Moderna 0.1mg/0.5mL dose dispersion for inj MDV",
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
        product_name_matches="COVID-19 mRNA (nucleoside modified) Vaccine Moderna 0.1mg/0.5mL dose dispersion for inj MDV",
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
        product_name_matches="COVID-19 mRNA (nucleoside modified) Vaccine Moderna 0.1mg/0.5mL dose dispersion for inj MDV",
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
    **covid_test_date_X(
        name="covid_test",
        index_date="elig_date + 43 days",
        n=6,
        test_result="any",
        return_expectations = {
            "date": {"earliest": start_date,  "latest" : end_date},
            "rate": "exponential_increase",
        },
    ),
    
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
    **covid_test_date_X(
        name="positive_test",
        index_date="elig_date + 43 days",
        n=6,
        test_result="positive",
        return_expectations={
            "date": {"earliest": start_date,  "latest" : end_date},
            "rate": "exponential_increase",
        },
    ),

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
    
    **with_these_clinical_events_date_X(
        name="primary_care_covid_case",
        n=6,
        index_date="elig_date + 43 days",
        codelist=covid_primary_care_probable_combined,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    
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
    **with_these_clinical_events_date_X(
        name = "primary_care_suspected_covid",
        n = 6,
        index_date = "elig_date + 43 days",
        codelist = primary_care_suspected_covid_combined,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    
    # emergency attendance
    **emergency_attendance_date_X(
        name = "emergency",
        n = 2,
        index_date = "elig_date + 43 days",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    
    # unplanned hospital admission
    admitted_unplanned_0_date=patients.admitted_to_hospital(
        returning="date_admitted",
        on_or_before="elig_date + 42 days",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    **admitted_date_X(
        name = "admitted_unplanned",
        n = 6,
        index_name = "admitted_unplanned",
        index_date = "elig_date + 43 days",
        returning="date_admitted",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    discharged_unplanned_0_date=patients.admitted_to_hospital(
        returning="date_discharged",
        on_or_before="elig_date + 42 days",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ), 
    **admitted_date_X(
        name = "discharged_unplanned",
        n = 6,
        index_name = "admitted_unplanned",
        index_date = "elig_date + 43 days",
        returning="date_discharged",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    
    # unplanned infectious hospital admission
    admitted_unplanned_infectious_0_date=patients.admitted_to_hospital(
        returning="date_admitted",
        on_or_before="elig_date + 42 days",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    **admitted_date_X(
        name = "admitted_unplanned_infectious",
        n = 6,
        index_name = "admitted_unplanned_infectious",
        index_date = "elig_date + 43 days",
        returning="date_admitted",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    discharged_unplanned_infectious_0_date=patients.admitted_to_hospital(
        returning="date_discharged",
        on_or_before="elig_date + 42 days",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ), 
    **admitted_date_X(
        name = "discharged_unplanned_infectious",
        n = 6,
        index_name = "admitted_unplanned_infectious",
        index_date = "elig_date + 43 days",
        returning="date_discharged",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    
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
    **admitted_date_X(
        name = "covidadmitted",
        n = 6,
        index_name = "covidadmitted",
        index_date = "elig_date + 42 days",
        returning="date_admitted",
        with_these_diagnoses=covid_codes,
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.05,
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

)