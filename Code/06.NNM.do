/*==============================================================================
    Project:    Treatment effects: NNM
    Author:     Shaun Scholes
    Date:       August 2026
    Purpose:    Estimate treatment effects in Stata
    Input:      cattaneo2.dta; rhc.xlsx
    Output:     Estimates
	Location:   "C:/Causal_Inference/Treatment_effects/code/06.NNM"
==============================================================================*/


clear all
cd "C:/Causal_Inference/"

*File location.
global projdir "./Treatment_effects/"
global data "${projdir}/data/"
global outputs "${projdir}/outputs/"
use "${data}/cattaneo2.dta", clear


gen id=_n
teffects nnmatch (bweight mmarried mage fage medu prenatal1) (mbsmoke), ///
ematch(mmarried prenatal1) nneighbor(1)  generate(mid)
teffects nnmatch (bweight mmarried mage fage medu prenatal1) (mbsmoke), ///
ematch(mmarried prenatal1) nneighbor(1) atet
di e(k_nnmax) 
local nmatch = e(k_nnmax) 

*link matches.
preserve
keep id bweight
sort id
save "${outputs}/bweight.dta",replace
restore

forvalues i = 1(1)`nmatch' {
preserve
keep mid`i'
gen id = mid`i'
merge m:1 id using "${outputs}/bweight.dta", noreport
keep if _merge==3
rename bweight mbweight`i'
keep mid`i' mbweight`i'
sort mid`i' 
save "${outputs}/Temp`i'.dta",replace
restore
}

forvalues i = 1(1)`nmatch' {
sort mid`i'	
merge m:m mid`i' using "${outputs}/Temp`i'.dta", nogen noreport
}

egen mbweight = rowmean(mbweight1-mbweight16)
*ATE = (Tx - control).
generate diff=0
replace diff = (mbweight - bweight) if mbsmoke==0        
replace diff = (bweight - mbweight) if mbsmoke==1       
qui: summ diff
di "ATE=" r(mean)
qui: summ diff if mbsmoke==1
di "ATET=" r(mean)

*erase intermediate files.

erase "${outputs}/bweight.dta"
forvalues i = 1/`nmatch' {
erase "${outputs}/Temp`i'.dta"
}

di "Finished"


******************
*NNM: Example 5
******************

clear
import excel "${data}/rhc.xlsx", sheet("rhc") firstrow

gen ARF=0
gen CHF=0
gen cirr=0
gen colcan=0
gen Coma=0
gen COPD=0
gen lungcan=0
gen MOSF=0
gen sepsis=0
gen female=0
gen died=0
gen treatment=0
replace ARF=1 if cat1=="ARF"
replace CHF=1 if cat1=="CHF"
replace cirr=1 if cat1=="Cirrhosis"
replace colcan=1 if cat1=="Colon Cancer"
replace Coma=1 if cat1=="Coma"
replace COPD=1 if cat1=="COPD"
replace lungcan=1 if cat1=="Lung Cancer"
replace MOSF=1 if cat1=="MOSF w/Malignancy"
replace sepsis=1 if cat1=="MOSF w/Sepsis"
replace female=1 if sex=="Female"
replace died=1 if death=="Yes"
replace treatment=1 if swang1=="RHC"
generate aps = aps1
keep ARF CHF cirr colcan Coma lungcan MOSF sepsis age female meanbp1 treatment died

*binary outcome
teffects nnmatch (died ARF CHF cirr colcan Coma lungcan MOSF ///
sepsis age female meanbp1) (treatment)

*similar set up but using IPW
teffects ipw (died) (treatment ARF CHF cirr colcan Coma lungcan ///
MOSF sepsis age female meanbp1, logit), ate


****************
**NNM: Example 6
****************

clear
import excel "${data}/rhc.xlsx", sheet("rhc") firstrow
gen ARF=0
gen CHF=0
gen cirr=0
gen colcan=0
gen Coma=0
gen COPD=0
gen lungcan=0
gen MOSF=0
gen sepsis=0
gen female=0
gen died=0
gen treatment=0
replace ARF=1 if cat1=="ARF"
replace CHF=1 if cat1=="CHF"
replace cirr=1 if cat1=="Cirrhosis"
replace colcan=1 if cat1=="Colon Cancer"
replace Coma=1 if cat1=="Coma"
replace COPD=1 if cat1=="COPD"
replace lungcan=1 if cat1=="Lung Cancer"
replace MOSF=1 if cat1=="MOSF w/Malignancy"
replace sepsis=1 if cat1=="MOSF w/Sepsis"
replace female=1 if sex=="Female"
replace died=1 if death=="Yes"
replace treatment=1 if swang1=="RHC"
generate aps = aps1
keep ARF CHF cirr colcan Coma lungcan MOSF sepsis age female meanbp1 treatment died
renvars, lower

set seed 68995
gen id=_n
*ATE
teffects nnmatch (died arf chf cirr colcan coma lungcan mosf ///
    sepsis age female meanbp1) (treatment), ///
gen(match) nneighbor(1) metric(mahalanobis) ate nolog

*ATET
teffects nnmatch (died arf chf cirr colcan coma lungcan mosf ///
sepsis age female meanbp1) (treatment), ///
nneighbor(1) metric(mahalanobis) atet nolog

sort id
generate age_of_match = age[match1]   /* 2184 matches  */
generate died_of_match = died[match1]  /* 2184 matches  */

*ATE
*IF everyone was on treatment
generate outcome_Tx=0
replace outcome_Tx=died if treatment==1
replace outcome_Tx=died_of_match if treatment==0
summ outcome_Tx
local POMt = r(mean)

*IF everyone was on control
generate outcome_control=0
replace outcome_control=died if treatment==0
replace outcome_control=died_of_match if treatment==1
summ outcome_control
local POMc = r(mean)
display (`POMt' - `POMc') 

*Calculate the ATT: observed Tx compared to matched controls
summ died if treatment==1       /* outcome for Tx */
local POMt = r(mean)
summ died_of_match if treatment==1    /* outcome for matched controls */
local POMc = r(mean)
display (`POMt' - `POMc')    

*=========================.
*FINISHED *********
*=========================.











