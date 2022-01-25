from cohortextractor import (
    StudyDefinition, 
    patients
)

# Import codelists.py script
from codelists import *

import pandas as pd

# import the vairables for deriving JCVI groups
from grouping_variables import (
    jcvi_variables, 
    study_parameters
)

# define variables explicitly from study_parameters
max_comparisons=study_parameters["max_comparisons"]
start_date=study_parameters["start_date"] # start of phase 1
end_date=study_parameters["end_date"] # latest date of data

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

# # import recurring event functions
# from recurrent_event_funs import *

# function for creating the comparison start dates for each elig_group
# avg_start_dates
avg_start_dates = pd.read_csv(
    filepath_or_buffer=f"./output/lib/avg_start_dates.csv",
    dtype=str
    )
# avg_start_dates for 2 ... K
def avg_start_k_date(K):
    
    def var_signature(k):
        return {
             f"avg_start_{k}_date": patients.categorised_as(
                    { avg_start_dates[f"avg_start_{k}_date"][i] : avg_start_dates['condition'][i] for i in avg_start_dates.index },
                    return_expectations={
                        "category":{
                            "ratios": { avg_start_dates[f"avg_start_{k}_date"][i] : 1/len(avg_start_dates.index) for i in avg_start_dates.index }
                        },
                    },
            ),
        }
    variables=dict()
    for k in range(1, K+1):
        variables.update(var_signature(k))
    return variables


# count tests in each of the K periods
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
            f"avg_start_{i}_date + {(i-1)*28 + 1} days", 
            f"avg_start_{i}_date + {i*28} days", 
            test_result, 
            return_expectations))
    return variables

# date of first test in each of the K periods
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
            f"avg_start_{i}_date + {(i-1)*28 + 1} days", 
            f"avg_start_{i}_date + {i*28} days", 
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

    # comparison start dates averaged over regions
    **avg_start_k_date(max_comparisons),

    # number of covid tests in each comparison period
    **covid_test_k_n(
        K=max_comparisons,
        test_result="any",
        return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
    ),
    # number of positive covid tests in each comparison period
    **covid_test_k_n(
        K=max_comparisons,
        test_result="positive",
        return_expectations={"int" : {"distribution": "poisson", "mean": 0.5}, "incidence" : 0.6}
    ),
    # first covid test in each comparison period
    **covid_test_k_date(
        K=max_comparisons,
        test_result="any",
        return_expectations={"date": {"earliest": start_date, "latest": end_date}}
    )

)