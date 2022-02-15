from cohortextractor import (
    StudyDefinition, 
    patients
)

# Import codelists.py script
from codelists import *

import pandas as pd

# import the vairables for deriving JCVI groups
from grouping_variables import (
    study_parameters
)

# define variables explicitly from study_parameters
max_comparisons=study_parameters["max_comparisons"]
start_date=study_parameters["start_date"] # start of phase 1
end_date=study_parameters["end_date"] # latest date of data

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

# import recurring event functions
from recurrent_event_funs import *

###
study=StudyDefinition(

    default_expectations={
        "date": {"earliest": start_date, "latest": end_date},
        "rate": "uniform",
        "incidence": 0.8,
    },  

    population=patients.all(),

    elig_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning='elig_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),
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
    # end date fo first comparison (start +28 days for vax, +56 days for unvax)
    end_1_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning='end_1_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),

    # age on start date
    age=patients.age_as_of("start_1_date"),

    ### covid tests as covariates
    # during unvaccinated time (from when tests widely availabe to elig_date)
    test_hist_1_n=patients.with_test_result_in_sgss(
                pathogen="SARS-CoV-2",
                test_result="any",
                between=["2020-05-18", "min_elig_date - 1 day"], # day before 1st vaccine eligibility date
                restrict_to_earliest_specimen_date=False,
                returning="number_of_matches_in_period",
                return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	        ),
    test_hist_2_n=patients.with_test_result_in_sgss(
                pathogen="SARS-CoV-2",
                test_result="any",
                between=["min_elig_date", "elig_date"],
                restrict_to_earliest_specimen_date=False,
                returning="number_of_matches_in_period",
                return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	        ),
    # during first dose time (elig date + 1 to elig_date + 6 weeks)
    test_hist_3_n=patients.with_test_result_in_sgss(
                pathogen="SARS-CoV-2",
                test_result="any",
                between=["elig_date + 1 day", "elig_date + 42 days"],
                restrict_to_earliest_specimen_date=False,
                returning="number_of_matches_in_period",
                return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	        ),
    
    ### covid tests as outcomes
    # number of covid tests in each 28-day period
    # **covid_test_k_n(
    #     K=max_comparisons,
    #     test_result="any",
    #     return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
    # ),
    # # number of positive covid tests in each 28-day period
    # **covid_test_k_n(
    #     K=max_comparisons,
    #     test_result="positive",
    #     return_expectations={"int" : {"distribution": "poisson", "mean": 0.5}, "incidence" : 0.6}
    # ),
    # first covid test in each 28-day period
    **covid_test_k_date(
        K=max_comparisons,
        test_result="any",
        return_expectations={"date": {"earliest": start_date, "latest": end_date}}
    ),

    #### clinical covariates
    # indicator for pregnancy at start of follow-up for comparison k
    **preg_36wks_k_date(max_comparisons),
    **pregdel_pre_k_date(max_comparisons),
    **preg_k(max_comparisons),

    bmi_0=patients.most_recent_bmi(
                on_or_before=f"start_1_date",
                minimum_age_at_measurement=16,
                include_measurement_date=True,
                date_format="YYYY-MM-DD"),
    **most_recent_bmi_k(max_comparisons),
    bmi_0_stage=patients.with_these_clinical_events(
                bmi_stage_primis,
                on_or_before=f"start_1_date",
                returning="category",
                include_date_of_match=True,
                date_format="YYYY-MM-DD",
                find_last_match_in_period=True
            ),
    **most_recent_bmi_stage_k(max_comparisons),

    # date of most recent shielding recording before each comparison start date
    shielded_0_date=patients.with_these_clinical_events(
        shield_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="start_1_date",
        find_last_match_in_period=True,
    ),
    **with_these_clinical_events_date_k(
        K=max_comparisons,
        name="shielded",
        codelist=shield_primis
    ),
    nonshielded_0_date=patients.with_these_clinical_events(
        shield_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="start_1_date",
        find_last_match_in_period=True,
    ),
    **with_these_clinical_events_date_k(
        K=max_comparisons,
        name="nonshielded",
        codelist=nonshield_primis
    ),
)