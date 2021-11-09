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
)