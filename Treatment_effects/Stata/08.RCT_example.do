/*==============================================================================
    Project:    Treatment effects: RCT_example
    Author:     Shaun Scholes
    Date:       August 2026
    Purpose:    Estimate treatment effects in Stata
    Input:      simulated
    Output:     Estimates
	Location:   "C:/Causal_Inference/Treatment_effects/code/08.RCT_example"
==============================================================================*/


**************************************************
*Three treatment groups, randomly assigned.
*compute means, differences in means, and t-test.
**************************************************

clear all


set obs 1000
set seed 66804457
gen id=_n
gen a = uniform()
gen Txgroup=0
replace Txgroup=1 if inrange(a,0.01,0.333)
replace Txgroup=2 if inrange(a,0.334,0.667)
replace Txgroup=3 if inrange(a,0.667,1)
label define Txgroup 1 "online" 2 "face-to-face" 3 "hybrid"
label values Txgroup Txgroup
tab1 Txgroup
gen b = uniform()
generate outcome = round(b*100,2)
bysort Txgroup: summ outcome
di (52.6875 - 41.73333)
bysort Txgroup: ci means outcome
ttest outcome if inlist(Txgroup,1,2), by(Txgroup)
keep id Txgroup outcome

bysort Txgroup: summ outcome
bysort Txgroup: ci means outcome
ttest outcome if inlist(Txgroup,1,2), by(Txgroup)

*=========================.
*FINISHED *********
*=========================.

