/*==============================================================================
    Project:    Treatment effects: PSmatching
    Author:     Shaun Scholes
    Date:       August 2026
    Purpose:    Estimate treatment effects in Stata
    Input:      cattaneo2.dta; psm.dta
    Output:     Estimates
	Location:   "C:/Causal_Inference/Treatment_effects/code/07.PS_matching"
==============================================================================*/


clear all
cd "C:/Causal_Inference/"

*File location.
global projdir "./Treatment_effects/"
global data "${projdir}/data/"
global outputs "${projdir}/outputs/"
set seed 67952
set obs 1000
gen id=_n
gen treat = uniform() < .4
gen outcome= runiformint(1,10)
gen x1= runiformint(1,10)
gen x2= runiformint(1,10)
teffects psmatch (outcome) (treat x1 x2,logit),gen(match)
teffects psmatch (outcome) (treat x1 x2,logit),atet
teffects psmatch (outcome) (treat x1 x2,logit),ate

local nmatch = e(k_nnmax) 

preserve
keep id outcome
sort id
save "${outputs}/outcome.dta",replace
restore

*max number of matches
local nmatch = e(k_nnmax)

forvalues i = 1(1)`nmatch' {
preserve
keep match`i'
gen id = match`i'
merge m:1 id using "${outputs}/outcome.dta", noreport
keep if _merge==3
rename outcome m_outcome`i'
keep match`i' m_outcome`i'
sort match`i' 
save "${outputs}/Temp`i'.dta",replace
restore
}

forvalues i = 1(1)`nmatch' {
sort match`i'	
merge m:m match`i' using "${outputs}/Temp`i'.dta", nogen noreport
}

egen m_outcome = rowmean(m_outcome1-m_outcome`nmatch')

*ATE = (Tx - control).
generate diff=0
replace diff = (m_outcome - outcome) if treat==0        
replace diff = (outcome - m_outcome) if treat==1       
qui: summ diff
di "ATE=" r(mean)
qui: summ diff if treat==1
di "ATET=" r(mean)


*erase intermediate files.
erase "${outputs}/outcome.dta"
forvalues i = 1/11 {
erase "${outputs}/Temp`i'.dta"
}


********************
*PSM: Example 
********************

clear
use "${data}/psm.dta", clear
 
rename t treat
*tab treat
*Unadjusted difference in means.
*ttest y, by(treat)
teffects psmatch (y) (treat x1 x2), atet
teffects psmatch (y) (treat x1 x2), gen(match) ate

*tebalance summarize

*We can add some predictions to see the propensity score (ps), 
*potential outcome (po), and treatment effect (te) for each observation.
predict ps0 ps1, ps
predict y0 y1, po
predict te
*Examine the first matched pair:
list treat y y0 y1 te if _n==1            /* control */
list treat y y0 y1 te if _n==467          /* matched to this Tx */
di (2.231719 - (-1.79457))                /* te for _n==1 (y1 - y0))*/
summarize te                              /* ATE */
summarize te if treat==1                  /* ATET */                 

*************************************************
*PSM: Example  (compare teffects and psmatch2)
*************************************************

clear
use "${data}/psm.dta", clear
rename t treat
teffects psmatch (y) (treat x1 x2), ate
teffects psmatch (y) (treat x1 x2), atet
*psmatch2 equivalent to teffects.
psmatch2 treat x1 x2, out(y) ate logit neighbor(1) 

*psmatch2 stores the estimate of the treatment effect on the treated in r(att), 
*this allows bootstrapping of the standard error of the 
*estimate (although it is unclear whether the
*bootstrap is valid in this context).  

bs "psmatch2 treat x1 x2, out(y)" "r(att)"


***************************
*2.5.3 Kernel matching 
***************************

clear
use "${data}/psm.dta", clear


rename t treat
psmatch2 treat x1 x2, kernel out(y) bwidth(0.5) common
bs "psmatch2 treat x1 x2, kernel out(y) bwidth(0.5) common" "r(att)", reps(100)
pstest x1 x2, t(treat) graph both


**PSM using kernel method**
*Our analysis used PSM, a technique that stimulates an experimental 
*setting in an observational data set and creates a treatment group and a 
*control group from the sample. One advantage of using PSM over regression approaches 
*is that it controls more effectively for the effects of observed confounders, 
*and hence while results remain observational, 
*bias attributable to confounding can be minimalised significantly. 
*We used PSM to estimate the average treatment effect for the treated (ATT), 
*which is the difference between the average mental health/wellbeing outcomes of 
*participants who had caring responsibilities (carers) 
*and the average outcomes for the same group under the hypothetical 
*scenario that they did not have any caring responsibilities (non-carers).
*Specifically, the PSM was performed on an unweighted data, 
*with the kernel matching method with 0.05 bandwidths to perform the matching. 
*Kernel matching uses weighted averages of all individuals in the control group to 
*create the counterfactual outcome, and matches participants in the treatment group to 
*those in the control groups based on the distance of their propensity score. 
*Higher weight is given to the matches whose propensity scores are closer to each other 
*and lower weight to those whose propensity scores are distal from each other. 
*A common support condition was imposed to ensure the quality of the matches. 
*95% confidence intervals were computed using bootstrapping 
*techniques with 100 replications. Missing values were handled with list-wise deletion. 
*High quality of matching was achieved; all analyses show Rubin's B<25%, 
*Rubin's R of 0.5-2, and a percentage bias of <10% for each covariate 
*This suggests that the unobservable heterogeneity reduced significantly after matching.


*******************************************
***psmatch2 
********************************************

*outcome: depress
*treatment: engagement
*matching variables: nssec wealth edqual worktime

clear all

set seed 67952
set obs 1000
gen id=_n
gen sex = uniform() < .5
gen depress = uniform() < .18
gen engage = uniform() < .18

gen nssec= runiformint(1,3)
gen wealth= runiformint(1,5)
gen educ= runiformint(1,3)
gen worktime= uniform() < .25

replace nssec=. if inrange(id,225,230)
replace wealth=. if inrange(id,230,239)
replace educ=. if inrange(id,330,339)
replace worktime=. if inrange(id,400,402)


*psmatch2 xvar, pscore(pvar) outcome(yvar) caliper(.001) noreplace neighbor(1)
psmatch2 engage sex nssec educ worktime, outcome(depress) caliper(.001) noreplace neighbor(1)
pstest educ sex nssec worktime, both treated(_treated) notable nodist label graph
summ _pscore   /* calculated by psmatch2: did not provide ourselves */

*187 matches (engage==1 plus matched engage==0)
*_id: In the case of one-to-one and NNM, a new identifier created for all observations.
summ _id   

*_nk In the case of one-to-one and NNM, 
*for every Tx observation, it stores the observation number of the k-th matched control observation. 
summ _n1 /* 187 (from 1 to 3709) */

*create pair id. 
gen id_pair = _id if _treated==0
replace id_pair = _n1 if _treated==1
bysort id_pair: egen paircount = count(id_pair)
drop if paircount !=2
count
label variable id_pair "(D) Pair identifier"
order id_pair id engage depress 

*depression for those with engage=0.
*depression for those with engage=1.

tab _treated

*mcc is used with matched case –control data.
*requires 1:1 matching (clogit for n:1)
*The data must be in the 1-observation-per-group format; that is, 
*the matched case and control must appear in 1 observation
*(the same format as required by the mcc command;

keep _treated depress id_pair
reshape wide depress, i(id_pair) j(_treated)

rename depress0 infrequent
rename depress1 frequent
tabulate frequent infrequent

*analysis of 187 matched pairs.

mcc frequent infrequent

di "Finished"


*=========================.
*FINISHED *********
*=========================.













