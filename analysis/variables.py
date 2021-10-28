
from datetime import date

from cohortextractor import (
    patients, 
)

# Import codelists.py script
from codelists import *

# import json module
import json

## import study dates
# change this in design.R if necessary
with open("./analysis/lib/dates.json") as f:
  studydates = json.load(f)

# define variables explicitly
ref_age_1=studydates["ref_age_1"] # reference date for calculating age for phase 1 groups
ref_age_2=studydates["ref_age_1"] # reference date for calculating age for phase 2 groups
ref_cev=studydates["ref_cev"] # reference date for calculating clinically extremely vulnerable group
ref_ar=studydates["ref_ar"] #reference date for caluclating at risk group
start_date=studydates["start_date"] # start of phase 1
end_date=studydates["end_date"] # end of followup
pandemic_start="2020-01-01"

## function to add days to a string date
from datetime import datetime, timedelta
def days(datestring, days):
  
  dt = datetime.strptime(datestring, "%Y-%m-%d").date()
  dt_add = dt + timedelta(days)
  datestring_add = datetime.strftime(dt_add, "%Y-%m-%d")

  return datestring_add

# Notes:
# for inequalities in the study definition, an extra expression is added to align with the comparison definitions in https://github.com/opensafely/covid19-vaccine-coverage-tpp-emis/blob/master/analysis/comparisons.py
# variables that define JCVI group membership MUST NOT be dependent on elig_date (index_date), this is for selecting the population based on registration dates and for deriving descriptive covariates
# JCVI groups are derived using ref_age_1, ref_age_2, ref_cev and ref_ar

jcvi_variables = dict(
  # age on phase 1 reference date
    age_1=patients.age_as_of(
        ref_age_1,
        return_expectations={
            "int": {"distribution": "population_ages"},
            "rate": "universal",
        },
    ),

    # age on phase 2 reference date
    age_2=patients.age_as_of(
        ref_age_2,
        return_expectations={
            "int": {"distribution": "population_ages"},
            "rate": "universal",
        },
    ),

    # patient sex
    sex=patients.sex(
        return_expectations={
        "rate": "universal",
        "category": {"ratios": {"M": 0.49, "F": 0.51}},
        "incidence": 1,
        }
    ),

    jcvi_group=patients.categorised_as(
        {
            "00": "DEFAULT",
            "01": "longres_group",
            "02": "age_1 >=80",
            "03": "age_1 >=75",
            "04": "age_1 >=70 OR (cev_group AND age_1 >=16 AND NOT preg_group)",
            "05": "age_1 >=65",
            "06": "atrisk_group AND age_1 >=16",
            "07": "age_1 >=60",
            "08": "age_1 >=55",
            "09": "age_1 >=50",
            "10": "age_2 >=40",
            "11": "age_2 >=30",
            "12": "age_2 >=18",
        },
        return_expectations={
            "rate": "universal",
            "incidence": 1,
            "category":{
                "ratios": {
                    "00": 1/13, "01": 1/13, "02": 1/13, "03": 1/13, "04": 1/13, "05": 1/13, "06": 1/13, "07": 1/13, "08":1/13, "09":1/13, "10":1/13, "11":1/13, "12":1/13}}
        },

    #### Pregnancy or Delivery codes recorded (for deriving JCVI group)
    # # date of last pregnancy code in 36 weeks before ref_cev
    preg_group=patients.satisfying(
        """
        (preg_36wks_date AND sex = 'F' AND age_1 < 50) AND
        (pregdel_pre_date <= preg_36wks_date OR NOT pregdel_pre_date)
        """,
        preg_36wks_date=patients.with_these_clinical_events(
            preg,
            returning="date",
            find_last_match_in_period=True,
            between=[days(ref_cev, -252), days(ref_cev, -1)],
            date_format="YYYY-MM-DD",
        ),
        # date of last delivery code recorded in 36 weeks before elig_date
        pregdel_pre_date=patients.with_these_clinical_events(
            pregdel,
            returning="date",
            find_last_match_in_period=True,
            between=[days(ref_cev, -252), days(ref_cev, -1)],
            date_format="YYYY-MM-DD",
        ),
    ),

    #### clinically extremely vulnerable group variables
    cev_group=patients.satisfying(
        "severely_clinically_vulnerable AND NOT less_vulnerable",

        # SHIELDED GROUP - first flag all patients with "high risk" codes
        severely_clinically_vulnerable=patients.with_these_clinical_events(
            shield,
            returning="binary_flag",
            on_or_before=days(ref_cev, -1),
            find_last_match_in_period=True,
        ),

        # find date at which the high risk code was added
        severely_clinically_vulnerable_date=patients.date_of(
            "severely_clinically_vulnerable",
            date_format="YYYY-MM-DD",
        ),

        # NOT SHIELDED GROUP (medium and low risk) - only flag if later than 'shielded'
        less_vulnerable=patients.with_these_clinical_events(
            nonshield,
            between=["severely_clinically_vulnerable_date + 1 day", days(ref_cev, -1)],
        ),
        return_expectations={"incidence": 0.01},
    ),

    #### at-risk group variables
    # Asthma Diagnosis code
    astdx=patients.with_these_clinical_events(
        ast,
        returning="binary_flag",
        on_or_before=days(ref_ar, -1),
        return_expectations={"incidence": 0.05},
    ),
            
    # asthma
    asthma_group=patients.satisfying(
        """
        astadm OR
        (astdx AND astrxm1 AND astrxm2 AND astrxm3)
        """,
        # day before date at which at risk group became eligible
        # Asthma Admission codes
        astadm=patients.with_these_clinical_events(
            astadm,
            returning="binary_flag",
            on_or_before=days(ref_ar, -1),
        ),
        # Asthma systemic steroid prescription code in month 1
        astrxm1=patients.with_these_medications(
            astrx,
            returning="binary_flag",
            between=[days(ref_ar, -31), days(ref_ar, -1)],
        ),
        # Asthma systemic steroid prescription code in month 2
        astrxm2=patients.with_these_medications(
            astrx,
            returning="binary_flag",
            between=[days(ref_ar, -61), days(ref_ar, -32)],
        ),
        # Asthma systemic steroid prescription code in month 3
        astrxm3=patients.with_these_medications(
            astrx,
            returning="binary_flag",
            between=[days(ref_ar, -91), days(ref_ar, -62)],
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Chronic Respiratory Disease other than asthma
    resp_group=patients.with_these_clinical_events(
        resp_cov,
        returning="binary_flag",
        on_or_before=days(ref_ar, -1),
        return_expectations={"incidence": 0.02},
    ),

    # Chronic Neurological Disease including Significant Learning Disorder
    cns_group=patients.with_these_clinical_events(
        cns_cov,
        returning="binary_flag",
        on_or_before=days(ref_ar, -1),
        return_expectations={"incidence": 0.01},
    ),

    # diabetes
    diab_group=patients.satisfying(
        """
        (NOT dmres_date AND diab_date) OR
        (dmres_date < diab_date)
        """,
        diab_date=patients.with_these_clinical_events(
            diab,
            returning="date",
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
            date_format="YYYY-MM-DD",
        ),
        dmres_date=patients.with_these_clinical_events(
            dmres,
            returning="date",
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
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
            sev_mental,
            returning="date",
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
            date_format="YYYY-MM-DD",
        ),
        # Remission codes relating to Severe Mental Illness
        smhres_date=patients.with_these_clinical_events(
            smhres,
            returning="date",
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
            date_format="YYYY-MM-DD",
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Chronic heart disease codes
    chd_group=patients.with_these_clinical_events(
        chd_cov,
        returning="binary_flag",
        on_or_before=days(ref_ar, -1),
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
            ckd15,
            returning="date",
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease codes-stages 3 - 5
        ckd35_date=patients.with_these_clinical_events(
            ckd35,
            returning="date",
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease diagnostic codes
        ckd=patients.with_these_clinical_events(
            ckd_cov,
            returning="binary_flag",
            on_or_before=days(ref_ar, -1),
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Chronic Liver disease codes
    cld_group=patients.with_these_clinical_events(
        cld,
        returning="binary_flag",
        on_or_before=days(ref_ar, -1),
        return_expectations={"incidence": 0.01},
    ),

    # immunosuppressed
    immuno_group=patients.satisfying(
        "immrx OR immdx", 
        # immunosuppression diagnosis codes
        immdx=patients.with_these_clinical_events(
            immdx_cov,
            returning="binary_flag",
            on_or_before=days(ref_ar, -1),
        ),
        # Immunosuppression medication codes
        immrx=patients.with_these_medications(
            immrx,
            returning="binary_flag",
            between=[days(ref_ar, -6*30), days(ref_ar, -1)],
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Asplenia or Dysfunction of the Spleen codes
    spln_group=patients.with_these_clinical_events(
        spln_cov,
        returning="binary_flag",
        on_or_before=days(ref_ar, -1),
        return_expectations={"incidence": 0.01},
    ),

    # Wider Learning Disability
    learndis_group=patients.with_these_clinical_events(
        learndis,
        returning="binary_flag",
        on_or_before=days(ref_ar, -1),
        return_expectations={"incidence": 0.01},
    ),

    # severe obesity
    sevobese_group=patients.satisfying(
        """
        (sev_obesity_date AND NOT bmi_date) OR
        (sev_obesity_date > bmi_date) OR
        bmi_value_temp >= 40
        """,
        bmi_stage_date=patients.with_these_clinical_events(
            bmi_stage,
            returning="date",
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
            date_format="YYYY-MM-DD",
        ),
        sev_obesity_date=patients.with_these_clinical_events(
            sev_obesity,
            returning="date",
            find_last_match_in_period=True,
            ignore_missing_values=True,
            between= ["bmi_stage_date", days(ref_ar, -1)],
            date_format="YYYY-MM-DD",
        ),
        bmi_date=patients.with_these_clinical_events(
            bmi,
            returning="date",
            ignore_missing_values=True,
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
            date_format="YYYY-MM-DD",
        ),
        bmi_value_temp=patients.with_these_clinical_events(
            bmi,
            returning="numeric_value",
            ignore_missing_values=True,
            find_last_match_in_period=True,
            on_or_before=days(ref_ar, -1),
            return_expectations={
                "float": {"distribution": "normal", "mean": 25, "stddev": 5},
            },
        ),
        return_expectations={"incidence": 0.01},
    ),

    # at risk group
    atrisk_group=patients.satisfying(
             """
             immuno_group OR
             ckd_group OR
             resp_group OR
             diab_group OR
             cld_group OR
             cns_group OR
             chd_group OR
             spln_group OR
             learndis_group OR
             sevment_group OR
             sevobese_group 
            """,
            return_expectations = {
            "incidence": 0.01,
            },
    ),

    # Patients in long-stay nursing and residential care
    longres_group=patients.with_these_clinical_events(
        longres,
        returning="binary_flag",
        on_or_before=days(start_date, -1),
        return_expectations={"incidence": 0.01},
    ),

    ### Will used the definition below, looks like it may be more thorough.
    # # CAREHOME STATUS
    # care_home_type=patients.care_home_status_as_of(
    #     "index_date - 1 day",
    #     categorised_as={
    #         "Carehome": """
    #           IsPotentialCareHome
    #           AND LocationDoesNotRequireNursing='Y'
    #           AND LocationRequiresNursing='N'
    #         """,
    #         "Nursinghome": """
    #           IsPotentialCareHome
    #           AND LocationDoesNotRequireNursing='N'
    #           AND LocationRequiresNursing='Y'
    #         """,
    #         "Mixed": "IsPotentialCareHome",
    #         "": "DEFAULT",  # use empty string
    #     },
    #     return_expectations={
    #         "category": {"ratios": {"Carehome": 0.05, "Nursinghome": 0.05, "Mixed": 0.05, "": 0.85, }, },
    #         "incidence": 1,
    #     },
    # ),

    # # simple care home flag
    # care_home_tpp=patients.satisfying(
    #     """care_home_type""",
    #     return_expectations={"incidence": 0.01},
    # ),
    
    # care_home_code=patients.with_these_clinical_events(
    #     carehome_primis_codes,
    #     on_or_before="index_date - 1 day",
    #     returning="binary_flag",
    #     return_expectations={"incidence": 0.01},
    # ),
    ),

    # vaccine eligibility dates
    elig_date=patients.categorised_as(
        {   ###
            "2020-12-08": "jcvi_group='01' OR jcvi_group='02' OR jcvi_group='03'",
            ###
            "2021-01-18": "jcvi_group='04'",
            ###
            "2021-02-15": "jcvi_group='05' OR jcvi_group='06'",
            ###
            "2021-02-22": "age_1 >= 64 AND age_1 < 65",
            "2021-03-01": "age_1 >= 60 AND age_1 < 64",
            ###
            "2021-03-08": "age_1 >= 56 AND age_1 < 60",
            "2021-03-09": "age_1 >= 55 AND age_1 < 56",
            ###
            "2021-03-19": "age_1 >= 50 AND age_1 < 55",
            ###
            "2021-04-13": "age_2 >= 45 AND age_1 < 50",
            "2021-04-26": "age_2 >= 44 AND age_1 < 45",
            "2021-04-27": "age_2 >= 42 AND age_1 < 44",
            "2021-04-30": "age_2 >= 40 AND age_1 < 42",
            ###
            "2021-05-13": "age_2 >= 38 AND age_2 < 40",
            "2021-05-19": "age_2 >= 36 AND age_2 < 38",
            "2021-05-21": "age_2 >= 34 AND age_2 < 36",
            "2021-05-25": "age_2 >= 32 AND age_2 < 34",
            "2021-05-26": "age_2 >= 30 AND age_2 < 32",
            ###
            "2021-06-08": "age_2 >= 25 AND age_2 < 30",
            "2021-06-15": "age_2 >= 23 AND age_2 < 25",
            "2021-06-16": "age_2 >= 21 AND age_2 < 23",
            "2021-06-18": "age_2 >= 18 AND age_2 < 21",
            "2100-12-31": "DEFAULT",
        },
        return_expectations={
            "category": {"ratios": 
            {
            "2020-12-08": 1/22,
            "2021-01-18": 1/22,
            "2021-02-15": 1/22,
            "2021-02-22": 1/22,
            "2021-03-01": 1/22,
            "2021-03-08": 1/22,
            "2021-03-09": 1/22,
            "2021-03-19": 1/22,
            "2021-04-13": 1/22,
            "2021-04-26": 1/22,
            "2021-04-27": 1/22,
            "2021-04-30": 1/22,
            "2021-05-13": 1/22,
            "2021-05-19": 1/22,
            "2021-05-21": 1/22,
            "2021-05-25": 1/22,
            "2021-05-26": 1/22,
            "2021-06-08": 1/22,
            "2021-06-15": 1/22,
            "2021-06-16": 1/22,
            "2021-06-18": 1/22,
            "2100-12-31": 1/22,
            }},
            "incidence": 1,
        },
    ),

)