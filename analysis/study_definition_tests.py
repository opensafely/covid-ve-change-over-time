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
# count tests in K comparison periods
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
    for i in range(1, K+1): 
        variables.update(var_signature(
            f"{test_result}_test_{i}_n", 
            f"start_1_date + {(i-1)*28 + 1} days", 
            f"end_1_date + {(i-1)*28} days", 
            test_result, 
            return_expectations))
    return variables

# date of first test in K comparison periods
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
    for i in range(1, K+1): 
        variables.update(var_signature(
            f"{test_result}_test_{i}_date", 
            f"start_1_date + {(i-1)*28 + 1} days", 
            f"end_1_date + {(i-1)*28} days", 
            test_result, 
            return_expectations))
    return variables

### pregnancy variables (derived in this script rather than study_definition.py as needed start_1_date)
# date of last pregnancy code in 36 weeks before start of comparison follow-up
def preg_36wks_k_date(K):

    def var_signature(name, lower, upper):
        return {
            name: patients.with_these_clinical_events(
                preg_primis,
                returning="date",
                find_last_match_in_period=True,
                between=[lower, upper],
                date_format="YYYY-MM-DD"),
        }
    variables = dict()
    for i in range(1, K+1): 
        variables.update(var_signature(
            name=f"preg_36wks_{i}_date",
            lower=f"start_1_date - {252 - (i-1)*28} days", 
            upper=f"start_1_date + {(i-1)*28} days"))
    return variables

# date of last delivery code recorded in 36 weeks before start of comparison follow-up
def pregdel_pre_k_date(K):

    def var_signature(name, lower, upper):
        return {
            name: patients.with_these_clinical_events(
            pregdel_primis,
            returning="date",
            find_last_match_in_period=True,
            between=[lower, upper],
            date_format="YYYY-MM-DD"),
        }
    variables = dict()
    for i in range(1, K+1): 
        variables.update(var_signature(
            name=f"pregdel_pre_{i}_date",
            lower=f"start_1_date - {252 - (i-1)*28} days", 
            upper=f"start_1_date + {(i-1)*28} days"))
    return variables

# pregnancy vairable that updates at start of each comparison
def preg_k(K):

    def var_signature(name, condition):
        return {
            name: patients.satisfying(condition),
        }
    variables = dict()
    for i in range(1, K+1): 
        variables.update(var_signature(
            name=f"preg_{i}", 
            condition=f"(preg_36wks_{i}_date) AND (pregdel_pre_{i}_date <= preg_36wks_{i}_date OR NOT pregdel_pre_{i}_date)"))
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
    ),

    # indicator for pregnancy at start of follow-up for comparison k
    **preg_36wks_k_date(K=max_comparisons),
    **pregdel_pre_k_date(K=max_comparisons),
    **preg_k(K=max_comparisons)

)