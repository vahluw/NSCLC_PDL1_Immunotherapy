**************************** 
/* Propensity score matching with all DAG factors --> but this time, kaplan-meier curves show days from treatment initiation until progression*/

/*wiresidence mnresidence inresidence varesidence prresidence dcresidence utresidence idresidence moresidence ctresidence nhresidence caresidence arresidence nvresidence deresidence mdresidence tnresidence alresidence njresidence paresidence nyresidence neresidence waresidence wvresidence waresidence azresidence laresidence orresidence okresidence txresidence coresidence iaresidence msresidence riresidence ohresidence scresidence garesidence miresidence ncresidence meresidence flresidence ilresidence nmresidence hiresidence ksresidence kyresidence maresidence */
****************************

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */

global indiv_covar_logistic "io_mono firstlinechemotherapy nonfirstlinechemotherapy otherfirstlinetherapy firstlinecombinationtherapy   antibrafdrug trkinhibitor metinhibitor carboplatinmonotherapy cisplatinmonotherapy    i.antialkdrug#i.alk i.antiegfrdrug#i.egfr i.antiros1drug#i.ros1 braf kras ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stage1 stage2 stage3 stage4 occult neversmoker  academicmedicalcenter kidney_bool liver_bool connective_tissue_bool interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused bevacizumabused threeormorechemotherapydrugs"

global psm "i.braf i.ecog0 i.ecog1 i.ecog2 i.ecog3 i.ecog4 i.nonsquamouscellcarcinoma  diagnosisyear ageatdiagnosis i.white i.asian i.black i.otherrace i.hispanicrace i.patientassistanceprogram i.othergovernmentalinsurance i.medicare i.selfpay i.medicaid i.commercialhealthplan i.noinsurance i.academicmedicalcenter i.kidney_bool i.liver_bool i.connective_tissue_bool i.interstitiallungdisease  i.bonemetastases i.brainmetastases   i.antiinfectiveusepriortotreatment i.glucocorticoidusepriortotreatmen pdl1 pdl1reported"


global indiv_covar_psm_prog_outcome "i.braf  i.ecog0 i.ecog1 i.ecog2 i.ecog3 i.ecog4 i.nonsquamouscellcarcinoma pdl1  diagnosisyear ageatdiagnosis i.white i.asian i.black i.otherrace i.hispanicrace i.patientassistanceprogram i.othergovernmentalinsurance i.medicare i.selfpay i.medicaid i.commercialhealthplan i.noinsurance  i.academicmedicalcenter i.bonemetastases i.brainmetastases  i.digestivesystemmetastases   i.antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen i.male"

global indiv_covar_psm_mort_outcome " i.ecog0 i.ecog1 i.ecog2 i.ecog3 i.ecog4  i.nonsquamouscellcarcinoma pdl1   diagnosisyear ageatdiagnosis i.white i.asian i.black i.otherrace i.hispanicrace i.male daysfromadvanceddiagnosistotreat i.patientassistanceprogram i.othergovernmentalinsurance i.medicare i.selfpay i.medicaid i.commercialhealthplan i.noinsurance i.neversmoker  i.academicmedicalcenter  i.diabetes i.bonemetastases i.brainmetastases i.othercnsmetastases i.digestivesystemmetastases i.adrenalmetastases i.unspecifiedmetastases   i.antiinfectiveusepriortotreatment i.glucocorticoidusepriortotreatmen i.braf  i.pdl1reported i.kidney_bool i.liver_bool i.connective_tissue_bool i.interstitiallungdisease"




global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
cd "${path}"
set scheme cleanplots




import delimited "all_data_365_10000.csv", clear 

tab progression_outcome
tab mortality_outcome
sum ageatdiagnosis
tab male
tab female
count if male==0 & female==0
tab white
tab black
tab asian
tab otherrace
tab hispanicethnicity
tab stage0
count if stagei==1 | stageia ==1 | stageia1 ==1 | stageia2 ==1 | stageia3 ==1 | stageib ==1 
gen stage1 = (stagei==1 | stageia ==1 | stageia1 ==1 | stageia2 ==1 | stageia3 ==1 | stageib ==1 )
count if stageii==1 | stageiia ==1 | stageiib ==1  
gen stage2 = (stageii==1 | stageiia ==1 | stageiib ==1)
count if stageiii ==1 | stageiiia ==1 | stageiiib | stageiiic ==1
gen stage3 = (stageiii ==1 | stageiiia ==1 | stageiiib | stageiiic ==1)
gen stage4 = (stageiv ==1 | stageiva ==1 | stageivb == 1)
count if stageiv ==1 | stageiva ==1 | stageivb == 1
count if occult==1
tab ecog0
tab ecog1
tab ecog2
tab ecog3
tab ecog4
count if ecog0==0 & ecog1==0 & ecog2==0 & ecog3 == 0 & ecog4==0
tab academicmedicalcenter
tab previoussmoker
tab neversmoker
tab squamouscellcarcinoma
tab nonsquamouscellcarcinoma
count if pdl1>=0.5
count if pdl1<0.5 & pdl1reported==1
count if pdl1reported==0
count if medicare==1 
count if medicaid==1
count if commercialhealthplan==1
count if noinsurance==1
count if selfpay==1
tab firstlinechemotherapy
tab firstlinecombinationtherapy
count if firstlinenivolumabmonotherapy ==1 | firstlinepembrolizumabmonotherap == 1 | firstlineatezolizumabmonotherapy == 1 | firstlinecemiplimabmonotherapy  == 1 | firstlinedurvalumabmonotherapy  == 1  | firstlineipilimumabnivolumab == 1
tab alk
tab egfr
tab ros1
tab egfr
tab kras
tab diabetes
gen kappa = 0.9
replace kappa = 0.7 if female==1
gen alpha = -0.302
replace alpha = -0.241 if female ==1 
gen extra_term = 1
replace extra_term = 1.012 if female==1
gen estimated_gfr = 142 * min(creatinine/kappa, 1)^(alpha) * max(creatinine/kappa, 1)^(-1.2) * 0.9938^(ageatdiagnosis) * extra_term

count if connectivetissuedisease ==1 | scleroderma == 1 | lupus== 1 | rheumatoidarthritis==1
count if chronickidneydisease == 1 |  priorkidneytransplant == 1 
count if interstitiallungdisease == 1
count if cirrhosis == 1 | hepatitis == 1 | priorlivertransplant  == 1

tab brainmetastases
tab bonemetastases
tab othercnsmetastases
tab digestivesystemmetastases
tab adrenalmetastases 
tab unspecifiedmetastases 
tab glucocorticoidusepriortotreatmen
tab antiinfectiveusepriortotreatment

gen connective_tissue_bool = ( connectivetissuedisease ==1 | scleroderma == 1 | lupus== 1 | rheumatoidarthritis==1)
gen kidney_bool = (chronickidneydisease == 1 |  priorkidneytransplant == 1 | (estimated_gfr < 60 & creatinine > 0))
gen liver_bool = (cirrhosis == 1 | hepatitis == 1 | priorlivertransplant  == 1 | ast >=109 | alt >=97 | bilirubin >= 2.0)
replace otherfirstlinetherapy = 1 if antirasdrug==1
replace daysfromadvanceddiagnosistotreat = log(daysfromadvanceddiagnosistotreat)
replace ageatdiagnosis = log(ageatdiagnosis)
replace diagnosisyear = 2022 - diagnosisyear

tab connective_tissue_bool
tab kidney_bool
tab liver_bool

gen time_limit = 365
gen outcome = "progression"

gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5

gen censor_time = progression_days 
replace censor_time = mortality_days if outcome == "mortality"
gen endpoint = progression_outcome
replace endpoint = mortality_outcome if outcome == "mortality"


tab progression_outcome
gen insured = 0
replace insured = 1 if patientassistanceprogram == 1 | othergovernmentalinsurance == 1 | medicare == 1 | medicaid == 1 | commercialhealthplan ==1


histogram pdl1, percent bin(9) xtitle("PD-L1 Intensity") color(ebblue) 
window manage forward results
graph export "PD_L1_distribution_whole_dataset_include_0_prog.png", replace

gen therapy_type = 1
replace therapy_type = 0 if firstlinechemotherapy == 1
replace therapy_type = 2 if firstlinecombinationtherapy == 1
replace therapy_type = 3 if nonfirstlinechemotherapy == 1
replace therapy_type = 4 if antialkdrug == 1
replace therapy_type = 5 if antiegfrdrug == 1
replace therapy_type = 6 if antibrafdrug == 1
replace therapy_type = 7 if antiros1drug == 1
replace therapy_type = 8  if otherfirstlinetherapy == 1 
replace therapy_type = 9 if trkinhibitor == 1
replace therapy_type = 10 if metinhibitor== 1
replace therapy_type = 11 if carboplatinmonotherapy == 1
replace therapy_type = 12 if cisplatinmonotherapy == 1

gen io_mono = (therapy_type==1)


replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit

putexcel set "progression_odds_all.xlsx", replace
logistic progression_outcome  ${indiv_covar_logistic}  pdl1reported  
putexcel (A1) = etable
putexcel clear

putexcel set "mortality_odds_all.xlsx", replace
logistic mortality_outcome   ${indiv_covar_logistic}   pdl1reported  
putexcel (A1) = etable
putexcel clear

putexcel set "progression_odds_only_pdl1.xlsx", replace
logistic progression_outcome  ${indiv_covar_logistic} io_mono#c.pdl1 nonfirstlinechemotherapy#c.pdl1 firstlinecombinationtherapy#c.pdl1 firstlinechemotherapy#c.pdl1 antialkdrug#c.pdl1 antiegfrdrug#c.pdl1 antiros1drug#c.pdl1 antibrafdrug#c.pdl1 trkinhibitor#c.pdl1 metinhibitor#c.pdl1 carboplatinmonotherapy#c.pdl1 cisplatinmonotherapy#c.pdl1 if pdl1reported==1
putexcel (A1) = etable
putexcel clear

putexcel set "mortality_odds_only_pdl1.xlsx", replace
logistic mortality_outcome  ${indiv_covar_logistic}   io_mono#c.pdl1 nonfirstlinechemotherapy#c.pdl1 firstlinecombinationtherapy#c.pdl1 firstlinechemotherapy#c.pdl1 antialkdrug#c.pdl1 antiegfrdrug#c.pdl1 antiros1drug#c.pdl1 antibrafdrug#c.pdl1 trkinhibitor#c.pdl1 metinhibitor#c.pdl1 carboplatinmonotherapy#c.pdl1 cisplatinmonotherapy#c.pdl1  if pdl1reported==1
putexcel (A1) = etable
putexcel clear



putexcel set "io_odds_only_pdl1_interaction.xlsx", replace
logistic therapy_type  ${covar_psm_logit_treatment} if pdl1reported==1
putexcel (A1) = etable
putexcel clear

drop if therapy_type > 2
drop if alk==1
drop if egfr==1
drop if ros1==1
drop if clinicalstudydrugused == 1 | bevacizumabused==1


stset censor_time, failure(endpoint)
stci, by(therapy_type) rmean
stcox therapy_type
sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation  (Days)") ytitle("Proportion at Risk") legend(order(0 "First-Line Chemotherapy" 1 "First-Line Chemotherapy" 2 "First-Line Combination Therapy")) 
window manage forward results
graph export "prog_survival_without_mutations_pre_match_no_combo_prog.png", replace
sts test therapy_type, logrank


tab therapy_type
drop if hispanicrace==1

capture teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm), vce(robust) control(0) osample(violator1)
drop if violator1==1
tab therapy_type

teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1reported==1, vce(robust) control(0)
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1reported==1, vce(robust) control(2)

/*
teffects ipw (progression_outcome) (therapy_type $psm) if pdl1reported==1, vce(robust) control(0)
teffects ipw (progression_outcome) (therapy_type $psm) if pdl1reported==1, vce(robust) control(2)
*/

tab therapy_type if pdl1>0
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1>0, vce(robust) control(0)
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1>0, vce(robust) control(2)

tab therapy_type if pdl1>=0.5
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1>=0.5, vce(robust) control(0)
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1>=0.5, vce(robust) control(2)


drop if digestivesystemmetastases==1 | ecog4==1
capture teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1>=0.01 & pdl1<0.5, vce(robust) control(0) 
tab therapy_type if pdl1>=0.01 & pdl1<0.5
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1>=0.01 & pdl1<0.5, vce(robust) control(0) 
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1>=0.01 & pdl1<0.5, vce(robust) control(2) 

drop if     1.therapy_type#1.selfpay != 0  |    1.therapy_type#1.brainmetastases != 0 
capture teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1reported==1 & pdl1==0, vce(robust) control(0) osample(violator3)
drop if violator3==1

tab therapy_type if pdl1reported==1 & pdl1==0
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1reported==1 & pdl1==0, vce(robust) control(0)
teffects ipwra (progression_outcome $indiv_covar_psm_prog_outcome, logit) (therapy_type $psm) if pdl1reported==1 & pdl1==0, vce(robust) control(2)
