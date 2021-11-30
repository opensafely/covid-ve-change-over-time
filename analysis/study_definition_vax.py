from cohortextractor import (
    StudyDefinition, 
    patients, 
    filter_codes_by_category
)

# Import codelists.py script
from codelists import *

import pandas as pd

# set seed so that dummy data can be reproduced
# np.random.seed(123456) 
# do this after create-study-definition-covs merged and store seed in study_parameters 

# import the vairables for deriving JCVI groups
from grouping_variables import (
    jcvi_variables, 
    start_date,
    end_date,
)

# regions
regions = pd.read_csv(
    filepath_or_buffer='./output/lib/regions.csv',
    dtype=str
)
ratio_regions = { regions['region'][i] : float(regions['ratio'][i]) for i in regions.index }

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
            on_or_before="elig_date + 42 days",
            returning="binary_flag",
        ),
    ),

    # Healthcare worker flag on vaccine record
    hscworker=patients.with_healthcare_worker_flag_on_covid_vaccine_record(
        returning="binary_flag",
        return_expectations={"incidence": 0.01},
        ),

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

    ethnicity_6_sus=patients.with_ethnicity_from_sus(
        returning="group_6",  
        use_most_frequent_code=True,
        return_expectations={
            "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
            "incidence": 0.8,
        },
    ),

    # IMD - quintile
    imd_0=patients.address_as_of(
            "elig_date + 42 days",
            returning="index_of_multiple_deprivation",
            round_to_nearest=100,
             return_expectations={
                         "category": {"ratios": {c: 1/320 for c in range(100,32100,100)}}
                         }
        ),

    # region - NHS England 9 regions
    region_0=patients.registered_practice_as_of(
        "elig_date + 42 days",
        returning="nuts1_region_name",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": ratio_regions,
            },
            "incidence": 0.99
        },
    ),

    ##########################
    ### CLINICAL VARIABLES ###
    ##########################

    # BMI
    bmi_0=patients.most_recent_bmi(
        on_or_before="elig_date + 42 days",
        minimum_age_at_measurement=16,
        # on_most_recent_day_of_measurement=True, # returning an error for some reason
        return_expectations={
            "float": {"distribution": "normal", "mean": 28, "stddev": 8},
            "incidence": 0.80,
        },
    ),

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

    ######################
    ### COVID VACCINES ###
    ######################

    # ## any covid vaccination, identified by target disease
    # covid_vax_disease_1_date=patients.with_tpp_vaccination_record(
    #     target_disease_matches="SARS-2 CORONAVIRUS",
    #     on_or_after=start_date,
    #     find_first_match_in_period=True,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     return_expectations={
    #         "date": {
    #             "earliest": start_date,  
    #             "latest": end_date,
    #         },
    #         "incidence": 0.5
    #     },
    # ),
    # covid_vax_disease_2_date=patients.with_tpp_vaccination_record(
    #     target_disease_matches="SARS-2 CORONAVIRUS",
    #     on_or_after="covid_vax_disease_1_date + 1 day",
    #     find_first_match_in_period=True,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     return_expectations={
    #         "date": {
    #             "earliest": start_date,  
    #             "latest": end_date,
    #         },
    #         "incidence": 0.5
    #     },
    # ),
    # covid_vax_disease_3_date=patients.with_tpp_vaccination_record(
    #     target_disease_matches="SARS-2 CORONAVIRUS",
    #     on_or_after="covid_vax_disease_2_date + 1 day",
    #     find_first_match_in_period=True,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #      return_expectations={
    #         "date": {
    #             "earliest": start_date,  
    #             "latest": end_date,
    #         },
    #         "incidence": 0.5
    #     },
    # ),

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
)