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

    # elig_date=patients.with_value_from_file(
    #     f_path='output/data/data_eligible_e.csv', 
    #     returning='elig_date', 
    #     returning_type='date',
    #     date_format='YYYY-MM-DD'
    #     ),
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
    # test_hist_2_n=patients.with_test_result_in_sgss(
    #             pathogen="SARS-CoV-2",
    #             test_result="any",
    #             between=["min_elig_date", "elig_date"],
    #             restrict_to_earliest_specimen_date=False,
    #             returning="number_of_matches_in_period",
    #             return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	#         ),
    # # during first dose time (elig date + 1 to elig_date + 6 weeks)
    # test_hist_3_n=patients.with_test_result_in_sgss(
    #             pathogen="SARS-CoV-2",
    #             test_result="any",
    #             between=["elig_date + 1 day", "elig_date + 42 days"],
    #             restrict_to_earliest_specimen_date=False,
    #             returning="number_of_matches_in_period",
    #             return_expectations={"int" : {"distribution": "poisson", "mean": 2}, "incidence" : 0.6}
	#         ),
    
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

    #### clinical covariates, defined at start of comparison period 1
    #### clinically extremely vulnerable group variables
    cev_group=patients.satisfying(

        "severely_clinically_vulnerable AND NOT less_vulnerable",

        # SHIELDED GROUP - first flag all patients with "high risk" codes
        severely_clinically_vulnerable=patients.with_these_clinical_events(
            shield_primis,
            returning="binary_flag",
            on_or_before="start_1_date",
            find_last_match_in_period=True,
        ),

        # find date at which the high risk code was added
        severely_clinically_vulnerable_date=patients.date_of(
            "severely_clinically_vulnerable",
            date_format="YYYY-MM-DD",
        ),

        # NOT SHIELDED GROUP (medium and low risk) - only flag if later than 'shielded'
        less_vulnerable=patients.with_these_clinical_events(
            nonshield_primis,
            between=["severely_clinically_vulnerable_date + 1 day", "start_1_date"],
        ),
        return_expectations={"incidence": 0.01},
    ),

    #### at-risk group variables
    # asthma
    asthma_group=patients.satisfying(
    """
      astadm OR
      (ast AND astrxm1 AND astrxm2 AND astrxm3)
      """,
    # Asthma Admission codes
    astadm=patients.with_these_clinical_events(
      astadm_primis,
      returning="binary_flag",
      on_or_before="start_1_date",
    ),
    # Asthma Diagnosis code
    ast = patients.with_these_clinical_events(
      ast_primis,
      returning="binary_flag",
      on_or_before="start_1_date",
    ),
    # Asthma systemic steroid prescription code in month 1
    astrxm1=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between=["start_1_date - 30 days", "start_1_date"],
    ),
    # Asthma systemic steroid prescription code in month 2
    astrxm2=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between=["start_1_date - 60 days", "start_1_date - 31 days"],
    ),
    # Asthma systemic steroid prescription code in month 3
    astrxm3=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between=["start_1_date - 90 days", "start_1_date - 61 days"],
    ),

  ),

    # Chronic Respiratory Disease other than asthma
    resp_group=patients.with_these_clinical_events(
        resp_primis,
        returning="binary_flag",
        on_or_before="start_1_date",
        return_expectations={"incidence": 0.02},
    ),

    # Chronic Neurological Disease including Significant Learning Disorder
    cns_group=patients.with_these_clinical_events(
        cns_primis,
        returning="binary_flag",
        on_or_before="start_1_date",
        return_expectations={"incidence": 0.01},
    ),

    # diabetes
    diab_group=patients.satisfying(
        """
        (NOT dmres_date AND diab_date) OR
        (dmres_date < diab_date)
        """,
        diab_date=patients.with_these_clinical_events(
            diab_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="start_1_date",
            date_format="YYYY-MM-DD",
        ),
        dmres_date=patients.with_these_clinical_events(
            dmres_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="start_1_date",
            date_format="YYYY-MM-DD",
        ),
        return_expectations={"incidence": 0.01},
    ),

    # severe mental illness codes
    sevment_group=patients.satisfying(
        """
        (NOT smhres_date AND sev_mental_date) OR
        smhres_date < sev_mental_date
        """,
        # Severe Mental Illness codes
        sev_mental_date=patients.with_these_clinical_events(
            sev_mental_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="start_1_date",
            date_format="YYYY-MM-DD",
        ),
        # Remission codes relating to Severe Mental Illness
        smhres_date=patients.with_these_clinical_events(
            smhres_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="start_1_date",
            date_format="YYYY-MM-DD",
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Chronic heart disease codes
    chd_group=patients.with_these_clinical_events(
        chd_primis,
        returning="binary_flag",
        on_or_before="start_1_date",
        return_expectations={"incidence": 0.01},
    ),

    # Chronic kidney disease diagnostic codes
    ckd_group=patients.satisfying(
        """
            ckd OR
            (ckd15_date AND 
            (ckd35_date >= ckd15_date) OR (ckd35_date AND NOT ckd15_date))
        """,
        # Chronic kidney disease codes - all stages
        ckd15_date=patients.with_these_clinical_events(
            ckd15_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="start_1_date",
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease codes-stages 3 - 5
        ckd35_date=patients.with_these_clinical_events(
            ckd35_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="start_1_date",
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease diagnostic codes
        ckd=patients.with_these_clinical_events(
            ckd_primis,
            returning="binary_flag",
            on_or_before="start_1_date",
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Chronic Liver disease codes
    cld_group=patients.with_these_clinical_events(
        cld_primis,
        returning="binary_flag",
        on_or_before="start_1_date",
        return_expectations={"incidence": 0.01},
    ),

    # immunosuppressed
    immuno_group=patients.satisfying(
        "immrx OR immdx", 
        # immunosuppression diagnosis codes
        immdx=patients.with_these_clinical_events(
            immdx_primis,
            returning="binary_flag",
            on_or_before="start_1_date",
        ),
        # Immunosuppression medication codes
        immrx=patients.with_these_medications(
            immrx_primis,
            returning="binary_flag",
            between=["start_1_date - 180 days", "start_1_date"],
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Asplenia or Dysfunction of the Spleen codes
    spln_group=patients.with_these_clinical_events(
        spln_primis,
        returning="binary_flag",
        on_or_before="start_1_date",
        return_expectations={"incidence": 0.01},
    ),

    # Wider Learning Disability
    learndis_group=patients.with_these_clinical_events(
        learndis_primis,
        returning="binary_flag",
        on_or_before="start_1_date",
        return_expectations={"incidence": 0.01},
    ),

    # indicator for pregnancy at start of follow-up for comparison k
    **preg_36wks_k_date(max_comparisons),
    **pregdel_pre_k_date(max_comparisons),
    **preg_k(max_comparisons),

    bmi_0=patients.most_recent_bmi(
                on_or_before=f"start_1_date",
                minimum_age_at_measurement=16,
                include_measurement_date=True,
                date_format="YYYY-MM-DD"),
    # **most_recent_bmi_k(max_comparisons),
    bmi_0_stage=patients.with_these_clinical_events(
                bmi_stage_primis,
                on_or_before=f"start_1_date",
                returning="category",
                include_date_of_match=True,
                date_format="YYYY-MM-DD",
                find_last_match_in_period=True
            ),
    # **most_recent_bmi_stage_k(max_comparisons),

    # date of most recent shielding recording before each comparison start date
    # shielded_0_date=patients.with_these_clinical_events(
    #     shield_primis,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="start_1_date",
    #     find_last_match_in_period=True,
    # ),
    # **with_these_clinical_events_date_k(
    #     K=max_comparisons,
    #     name="shielded",
    #     codelist=shield_primis
    # ),
    # nonshielded_0_date=patients.with_these_clinical_events(
    #     shield_primis,
    #     returning="date",
    #     date_format="YYYY-MM-DD",
    #     on_or_before="start_1_date",
    #     find_last_match_in_period=True,
    # ),
    # **with_these_clinical_events_date_k(
    #     K=max_comparisons,
    #     name="nonshielded",
    #     codelist=nonshield_primis
    # ),

)