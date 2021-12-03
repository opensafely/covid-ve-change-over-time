from cohortextractor import (
    patients,
)

# define recurrent event variables

# imd defined at each start_k_date 
# (was not working with categorised_as, so do the categorising in R script)
def imd_k(K):
    
    def var_signature(name, date):
        return {
            name: patients.address_as_of(
                    date,
                    returning="index_of_multiple_deprivation",
                    round_to_nearest=100,
                    return_expectations={
                        "category": {"ratios": {c: 1/320 for c in range(100,32100,100)}}
                        }
                    ),
        }
    variables = {}
    for k in range(1, K+1):
        variables.update(var_signature(f"imd_{k}", f"start_1_date + {(k-1)*28} days"))
    return variables

# region defined at each start_k_date
def region_k(K):
    
    def var_signature(name, date):
        return {
            name: patients.registered_practice_as_of(
                date=date,
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
        }
    variables=var_signature(f"region_1", "start_1_date")
    for k in range(2, K+1):
        variables.update(var_signature(f"region_{k}", f"start_1_date + {(k-1)*28} days"))
    return variables


# bmi
def most_recent_bmi_X(name, index_date, n, return_expectations):

    def var_signature(name, on_or_before, return_expectations):
        return {
            name: patients.most_recent_bmi(
                on_or_before=on_or_before,
                minimum_age_at_measurement=16,
                include_measurement_date=True,
                date_format="YYYY-MM-DD",
                return_expectations=return_expectations
            )
        }
    variables=var_signature(f"{name}_1", index_date, return_expectations)
    for i in range(2, n+1):
        variables.update(var_signature(f"{name}_{i}", f"{name}_{i-1}_date_measured - 1 day", return_expectations))
    return variables


# clinical events with codelist
def with_these_clinical_events_date_X(name, codelist, index_date, n, return_expectations):
    
    def var_signature(name, on_or_after, codelist, return_expectations):
        return {
            name: patients.with_these_clinical_events(
                    codelist,
                    returning="date",
                    on_or_after=on_or_after,
                    date_format="YYYY-MM-DD",
                    find_first_match_in_period=True,
                    return_expectations=return_expectations
	        ),
        }
    variables=var_signature(f"{name}_1_date", index_date, codelist, return_expectations)
    for i in range(2, n+1):
        variables.update(var_signature(f"{name}_{i}_date", f"{name}_{i-1}_date + 1 day", codelist, return_expectations))
    return variables

# # medications with codelist
# def with_these_medications_date_X(name, codelist, index_date, n, return_expectations):
    
#     def var_signature(name, on_or_after, codelist, return_expectations):
#         return {
#             name: patients.with_these_medications(
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

# # bmi
# def most_recent_bmi_X(name, index_date, n, return_expectations):

#     def var_signature(name, on_or_before, return_expectations):
#         return {
#             name: patients.most_recent_bmi(
#                 on_or_before=on_or_before,
#                 minimum_age_at_measurement=16,
#                 include_measurement_date=True,
#                 date_format="YYYY-MM-DD",
#                 return_expectations=return_expectations
#             )
#         }
#     variables=var_signature(f"{name}_1", index_date, return_expectations)
#     for i in range(2, n+1):
#         variables.update(var_signature(f"{name}_{i}", f"{name}_{i-1}_date_measured - 1 day", return_expectations))
#     return variables

# # covid test date
# def covid_test_date_X(name, index_date, n, test_result, return_expectations):
    
#     def var_signature(name, on_or_after, test_result, return_expectations):
#         return {
#             name: patients.with_test_result_in_sgss(
#                 pathogen="SARS-CoV-2",
#                 test_result=test_result,
#                 on_or_after=on_or_after,
#                 find_first_match_in_period=True,
#                 restrict_to_earliest_specimen_date=False,
#                 returning="date",
#                 date_format="YYYY-MM-DD",
#                 return_expectations=return_expectations
# 	        ),
#         }
#     variables = var_signature(f"{name}_1_date", index_date, test_result, return_expectations)
#     for i in range(2, n+1):
#         variables.update(var_signature(f"{name}_{i}_date", f"{name}_{i-1}_date + 1 day", test_result, return_expectations))
#     return variables

# # emergency attendence
# def emergency_attendance_date_X(name, index_date, n, return_expectations):
    
#     def var_signature(name, on_or_after, return_expectations):
#         return {
#             name: patients.attended_emergency_care(
#                     returning="date_arrived",
#                     on_or_after=on_or_after,
#                     find_first_match_in_period=True,
#                     date_format="YYYY-MM-DD",
#                     return_expectations=return_expectations
# 	        ),
#         }
#     variables = var_signature(f"{name}_1_date", index_date, return_expectations)
#     for i in range(2, n+1):
#         variables.update(var_signature(f"{name}_{i}_date", f"{name}_{i-1}_date + 1 day", return_expectations))
#     return variables

# # hospital admissions
# def admitted_date_X(
#     name, index_name, index_date, n, returning, 
#     with_these_diagnoses=None, 
#     with_admission_method=None, 
#     with_patient_classification=None, 
#     return_expectations=None
# ):
#     def var_signature(
#         name, on_or_after, returning, 
#         with_these_diagnoses, 
#         with_admission_method, 
#         with_patient_classification, 
#         return_expectations
#     ):
#         return {
#             name: patients.admitted_to_hospital(
#                     returning = returning,
#                     on_or_after = on_or_after,
#                     find_first_match_in_period = True,
#                     date_format = "YYYY-MM-DD",
#                     with_these_diagnoses = with_these_diagnoses,
#                     with_admission_method = with_admission_method,
#                     with_patient_classification = with_patient_classification,
#                     return_expectations = return_expectations
# 	        ),
#         }
#     variables = var_signature(
#         f"{name}_1_date", 
#         index_date, 
#         returning, 
#         with_these_diagnoses,
#         with_admission_method,
#         with_patient_classification,
#         return_expectations
#     )
#     for i in range(2, n+1):
#         variables.update(var_signature(
#             f"{name}_{i}_date", f"{index_name}_{i-1}_date + 1 day", returning, 
#             with_these_diagnoses,
#             with_admission_method,
#             with_patient_classification,
#             return_expectations
#         ))
#     return variables