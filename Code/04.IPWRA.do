/*==============================================================================
    Project:    Treatment effects: IPWRA
    Author:     Shaun Scholes
    Date:       August 2026
    Purpose:    Estimate treatment effects in Stata
    Input:      cattaneo2.dta
    Output:     Estimates
	Location:   "C:/Causal_Inference/Treatment_effects/code/04.IPWRA"
==============================================================================*/


clear all
cd "C:/Causal_Inference/"

*File location.
global projdir "./Treatment_effects/"
global data "${projdir}/data/"
global outputs "${projdir}/outputs/"
use "${data}/cattaneo2.dta", clear

qui: teffects ipwra (bweight fage mage mmarried prenatal1) ///
(mbsmoke fage mage mmarried fbaby,logit), pomeans aequations nolog
teffects ipwra (bweight fage mage mmarried prenatal1) ///
(mbsmoke fage mage mmarried fbaby,logit), ate aequations nolog

qui: logit mbsmoke fage mage mmarried fbaby
predict pscore
generate ipw=0
replace ipw=1/(1-pscore) if mbsmoke==0
replace ipw=1/(pscore) if mbsmoke==1
qui: reg bweight fage mage mmarried prenatal1 if mbsmoke==1 [pweight=ipw]
predict double pom_t
qui: reg bweight fage mage mmarried prenatal1 if mbsmoke==0 [pweight=ipw]
predict double pom_c
qui: summ pom_t
local pom_t = r(mean)
qui: summ pom_c
local pom_c = r(mean)
di "ATE="`pom_t' - `pom_c'


di "Finished"


*=========================.
*FINISHED *********
*=========================.




