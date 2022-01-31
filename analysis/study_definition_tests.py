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

# recurrent var functions
# count tests in K 28-day periods
def covid_test_k_n(K, test_result, return_expectations):
    
    def var_signature(name, lower, upper, test_result, return_expectations):
        return {
            name: patients.with_test_result_in_sgss(
                pathogen="SARS-CoV-2",
                test_result=test_result,
                between=[lower, upper],
                restrict_to_earliest_specimen_date=False,
                returning="number_of_matches_in_period",
                return_expectations=return_expectations
	        ),
        }
    variables = dict()
    for i in range(1, K+2): # last one won't be used for the vax arm
        variables.update(var_signature(
            f"{test_result}_test_{i}_n", 
            f"start_1_date + {(i-1)*28 + 1} days", 
            f"start_1_date + {i*28} days", 
            test_result, 
            return_expectations))
    return variables

# date of first test in K 28-day periods
def covid_test_k_date(K, test_result, return_expectations):
    
    def var_signature(name, lower, upper, test_result, return_expectations):
        return {
            name: patients.with_test_result_in_sgss(
                pathogen="SARS-CoV-2",
                test_result=test_result,
                between=[lower, upper],
                restrict_to_earliest_specimen_date=False,
                find_first_match_in_period=True,
                returning="date",
                date_format = "YYYY-MM-DD",
                return_expectations=return_expectations
	        ),
        }
    variables = dict()
    for i in range(1, K+2): # last one won't be used for the vax arm
        variables.update(var_signature(
            f"{test_result}_test_{i}_date", 
            f"start_1_date + {(i-1)*28 + 1} days", 
            f"start_1_date + {i*28} days", 
            test_result, 
            return_expectations))
    return variables

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

    start_1_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning='start_1_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),

    ### covid tests as covariates
    # during unvaccinated time (from when tests widely availabe to elig_date)
    covid_test_pre_elig_n = patients.with_test_result_in_sgss(
                pathogen="SARS-CoV-2",
                test_result="any",
                between=["2020-05-18", "elig_date"],
                restrict_to_earliest_specimen_date=False,
                returning="number_of_matches_in_period",
                return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	        ),
    # during first dose time (elig date + 1 to elig_date + 6 weeks)
    covid_test_post_elig_n = patients.with_test_result_in_sgss(
                pathogen="SARS-CoV-2",
                test_result="any",
                between=["elig_date + 1 day", "elig_date + 42 days"],
                restrict_to_earliest_specimen_date=False,
                returning="number_of_matches_in_period",
                return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	        ),
    
    ### covid tests as outcomes
    # number of covid tests in each 28-day period
    **covid_test_k_n(
        K=max_comparisons,
        test_result="any",
        return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
    ),
    # number of positive covid tests in each 28-day period
    **covid_test_k_n(
        K=max_comparisons,
        test_result="positive",
        return_expectations={"int" : {"distribution": "poisson", "mean": 0.5}, "incidence" : 0.6}
    ),
    # first covid test in each 28-day period
    **covid_test_k_date(
        K=max_comparisons,
        test_result="any",
        return_expectations={"date": {"earliest": start_date, "latest": end_date}}
    )

)