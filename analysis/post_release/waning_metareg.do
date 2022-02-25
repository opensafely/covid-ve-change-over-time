cd "/Users/eh1415/Documents/covid-change-ve-over-time"

capture log close
log using waning_metareg, replace

set more off
clear


import delimited "/Users/eh1415/Documents/covid-ve-change-over-time/release20220221/data_metareg.csv"

replace estimate="" if estimate=="NA"
destring estimate, replace

replace confhigh="" if confhigh=="NA"
destring confhigh, replace

replace conflow="" if conflow=="NA"
destring conflow, replace

gen loghr=log(estimate)
gen lnu=log(confhigh)
gen lnl=log(conflow)

gen se1=(lnu-loghr)/1.96
gen se2=(loghr-lnl)/1.96

corr se1 se2

gen seloghr=(se1+se2)/2
drop se1 se2 lnu lnl

gen stratum=0 if subgroup=="65+ years"
replace stratum=1 if subgroup=="16-64 years and clinically vulnerable"
replace stratum=2 if subgroup=="40-64 years"
replace stratum=3 if subgroup=="18-39 years"

rename outcome temp

gen outcome=1 if temp=="covidadmitted"
replace outcome=2 if temp=="coviddeath"
replace outcome=3 if temp=="postest"
replace outcome=4 if temp=="noncoviddeath"
replace outcome=5 if temp=="anytest"

label define outcome 1 "COVID-19 hospitalisation" 2 "COVID-19 death" 3 "Positive test" ///
 4 "Non-COVID death" 5  "Any test"
label values outcome outcome
drop temp

encode comparison, gen(vaccine)
tab vaccine
label list vaccine

replace k=k-1

sort outcome stratum vaccine
save "/Users/eh1415/Documents/covid-ve-change-over-time/release20220221/waning_metareg.dta", replace
/*
metareg loghr k if outcome==1 & stratum==0 & vaccine==3, wsse(seloghr)
local a=_b[k]
local b=_se[k]
local c=_b[_cons]
local d=_se[_cons]

di `a'
di `b'
di `c'
di `d'

tempname memhold
postfile `memhold' outcome stratum vaccine logrhr selogrhr loghr1 seloghr1 using results, replace

  forvalues i=1/5 {
  	forvalues v=1/3 {
		forvalues s=0/3 {
				di "A: " `i' `s' `v'

				count if outcome==`i' & stratum==`s' & vaccine==`v' &loghr<.
				if r(N)>2 {
				di "B: " `i' `s' `v'
				metareg loghr k if outcome==`i' & stratum==`s' & vaccine==`v', wsse(seloghr)
				local a=_b[k]
				local b=_se[k]
				local c=_b[_cons]
				local d=_se[_cons]
				post `memhold' (`i') (`s') (`v') (`a') (`b') (`c') (`d')
				}
		}
	}
  }

postclose `memhold'
di "`memhold'"
*/
use results, clear
sort outcome stratum vaccine
save results, replace

use waning_metareg, clear
merge m:1 outcome stratum vaccine using results
tab _merge

drop if _merge<3
drop _merge
drop estimate conflow confhigh loghr seloghr
sort outcome stratum vaccine k
by outcome stratum vaccine: keep if _n==_N
replace k=k+1
assert k==6

drop stratum k vaccine

order subgroup comparison outcome model
*export excel using metareg_results.xlsx, replace firstrow(variables)

gen rhr=exp(logrhr)
gen lci=exp(logrhr-1.96*selogrhr)
gen uci=exp(logrhr+1.96*selogrhr)

encode comparison, gen(vaccine)
tab vaccine
label list vaccine


drop loghr1 seloghr1 comparison
sort outcome subgroup vaccine

describe

reshape wide logrhr selogrhr rhr lci uci, i(subgroup outcome model) j(vaccine)

gen stratum=0 if subgroup=="65+ years"
replace stratum=1 if subgroup=="16-64 years and clinically vulnerable"
replace stratum=2 if subgroup=="40-64 years"
replace stratum=3 if subgroup=="18-39 years"

label define stratum 0 "65+_years" 1 "16-64_years_&_clinically_vulnerable" ///
2 "40-64_years" 3 "18-39_years"

label values stratum stratum

sort outcome stratum
compress

foreach var in "rhr" "lci" "uci" {
	format `var'1 `var'2 `var'3 %5.2f
	rename `var'1 `var'_PB
	rename `var'2 `var'_AZ
	rename `var'3 `var'_PB_vs_AZ
}

*foreach vacc in "PB" "AZ" "PV_vs_AZ"

describe
list outcome stratum rhr_PB lci_PB uci_PB rhr_AZ lci_AZ uci_AZ if outcome<4, noobs nodisp clean

list outcome stratum rhr_PB lci_PB uci_PB rhr_AZ lci_AZ uci_AZ if outcome>3, noobs nodisp clean

list outcome stratum rhr_PB_vs_AZ lci_PB_vs_AZ uci_PB_vs_AZ if outcome<4 & rhr_PB_vs_AZ<., noobs nodisp clean

list outcome stratum rhr_PB_vs_AZ lci_PB_vs_AZ uci_PB_vs_AZ if outcome>3 & rhr_PB_vs_AZ<., noobs nodisp clean



