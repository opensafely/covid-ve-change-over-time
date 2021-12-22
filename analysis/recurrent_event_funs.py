from cohortextractor import (
    patients,
)

# define recurrent event variables
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

# hospital admissions
def admitted_date_X(
    name, index_name, index_date, n, returning, 
    with_these_diagnoses=None, 
    with_admission_method=None, 
    with_patient_classification=None, 
    return_expectations=None
):
    def var_signature(
        name, on_or_after, returning, 
        with_these_diagnoses, 
        with_admission_method, 
        with_patient_classification, 
        return_expectations
    ):
        return {
            name: patients.admitted_to_hospital(
                    returning = returning,
                    on_or_after = on_or_after,
                    find_first_match_in_period = True,
                    date_format = "YYYY-MM-DD",
                    with_these_diagnoses = with_these_diagnoses,
                    with_admission_method = with_admission_method,
                    with_patient_classification = with_patient_classification,
                    return_expectations = return_expectations
	        ),
        }
    variables = var_signature(
        f"{name}_1_date", 
        index_date, 
        returning, 
        with_these_diagnoses,
        with_admission_method,
        with_patient_classification,
        return_expectations
    )
    for i in range(2, n+1):
        variables.update(var_signature(
            f"{name}_{i}_date", f"{index_name}_{i-1}_date + 1 day", returning, 
            with_these_diagnoses,
            with_admission_method,
            with_patient_classification,
            return_expectations
        ))
    return variables