from cohortextractor import (
    patients,
)

# Import codelists.py script
from codelists import *

# # clinical events with codelist
# def with_these_clinical_events_date_X(name, codelist, index_date, n, return_expectations):
    
#     def var_signature(name, on_or_after, codelist, return_expectations):
#         return {
#             name: patients.with_these_clinical_events(
#                     codelist,
#                     returning="date",
#                     on_or_after=on_or_after,
#                     date_format="YYYY-MM-DD",
#                     find_first_match_in_period=True,
#                     return_expectations=return_expectations
# 	        ),
#         }
#     variables=var_signature(f"{name}_1_date", index_date, codelist, return_expectations)
#     for i in range(2, n+1):
#         variables.update(var_signature(f"{name}_{i}_date", f"{name}_{i-1}_date + 1 day", codelist, return_expectations))
#     return variables

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

### pregnancy variables 
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

# date of:
# - most recent numeric BMI before start_1_date
# - latest BMI numeric recording during each comparison period
def most_recent_bmi_k(K):

    def var_signature(name, lower, upper):
        return {
            name: patients.most_recent_bmi(
                between=[lower,upper],
                minimum_age_at_measurement=16,
                include_measurement_date=True,
                date_format="YYYY-MM-DD"
            )
        }
    variables=dict()
    for i in range(1, K):
        variables.update(var_signature(
            name=f"bmi_{i}", 
            lower=f"start_1_date + {(i-1)*28 + 1} days", 
            upper=f"end_1_date + {(i-1)*28} days"))
    return variables

def most_recent_bmi_stage_k(K):

    def var_signature(name, lower, upper):
        return {
            name: patients.with_these_clinical_events(
                bmi_stage_primis,
                between=[lower,upper],
                returning="category",
                include_date_of_match=True,
                date_format="YYYY-MM-DD",
                find_last_match_in_period=True
            )
        }
    variables=dict()
    for i in range(1, K):
        variables.update(var_signature(
            name=f"bmi_{i}_stage", 
            lower=f"start_1_date + {(i-1)*28 + 1} days", 
            upper=f"end_1_date + {(i-1)*28} days"))
    return variables


# clinical events with codelist in each comparison period
def with_these_clinical_events_date_k(K, name, codelist):
    
    def var_signature(name, codelist, lower, upper):
        return {
            name: patients.with_these_clinical_events(
                    codelist,
                    returning="date",
                    between=[lower,upper],
                    date_format="YYYY-MM-DD",
                    find_last_match_in_period=True,
	        ),
        }
    variables=dict()
    for i in range(1, K):
        variables.update(var_signature(
            name=f"{name}_{i}_date", 
            codelist=codelist,
            lower=f"start_1_date + {(i-1)*28 + 1} days", 
            upper=f"end_1_date + {(i-1)*28} days"),
            )
    return variables
