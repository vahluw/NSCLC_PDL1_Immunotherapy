global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/clf_xgb/"
cd "${path}"
set scheme cleanplots

import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/clf_xgb/test_set_365_10000_1.csv", clear 
rocreg progression_outcome prog_preds

rocreg mortality_outcome mort_preds

import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/clf_xgb/test_set_365_10000_1.csv", clear 

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

replace daysfromadvanceddiagnosistotreat = log(daysfromadvanceddiagnosistotreat)
replace ageatdiagnosis = log(ageatdiagnosis)
replace diagnosisyear = 2022 - diagnosisyear

tab connective_tissue_bool
tab kidney_bool
tab liver_bool



gen threshold = 0.686
gen time_limit = 365


gen prog_pred = prog_preds
gen pdl1_over_threshold = (pdl1 >=0.5) 

replace diagnosisyear = 2024 - diagnosisyear
gen progressed_prediction = (prog_pred>=threshold)
gen hazard_prediction  = 0
replace hazard_prediction = 1 if progressed_prediction==0

sum progression_days
replace progression_days = time_limit if progression_days == 0 | progression_days >time_limit
replace mortality_days = time_limit if mortality_days == 0 | mortality_days >time_limit



capture stset progression_days, failure(progression_outcome)
stci, by(hazard_prediction) rmean
stcox hazard_prediction
sts graph, by(hazard_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "test_set_ml_mutations_xgb_prog_io.png", replace


drop if alk==1
drop if egfr==1
drop if ros1==1


capture stset progression_days if pdl1reported==1, failure(progression_outcome)
stci if  pdl1reported==1, by(pdl1_over_threshold) rmean
stcox pdl1_over_threshold if  pdl1reported==1
sts graph if  pdl1reported==1, by(pdl1_over_threshold) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by PD-L1 Status; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "PD-L1 <50%" 2 "PD-L1 >=50%"))
graph export "io_test_set_pdl1_prog_io.png", replace

capture stset progression_days, failure(progression_outcome)
stci, by(hazard_prediction) rmean
stcox hazard_prediction
sts graph, by(hazard_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "test_set_ml_no_mutations_xgb_prog_io.png", replace

capture stset progression_days if  pdl1reported==1, failure(progression_outcome)
stci if  pdl1reported==1, by(hazard_prediction) rmean 
stcox hazard_prediction if  pdl1reported==1
sts graph, by(hazard_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "io_only_test_set_xgb_prog_pdl1reported_io.png", replace


capture stset progression_days if  pdl1>=0.5, failure(progression_outcome)
stci if  pdl1>=0.5, by(hazard_prediction) rmean 
stcox hazard_prediction if  pdl1>=0.5
sts graph if  pdl1>=0.5, by(hazard_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "io_only_test_set_pdl1_over50_xgb_prog_io.png", replace

capture stset progression_days if  pdl1<0.5 & pdl1reported==1, failure(progression_outcome)
stci if  pdl1<0.5 & pdl1reported==1, by(hazard_prediction) rmean 
stcox hazard_prediction if  pdl1<0.5 & pdl1reported==1
sts graph if  pdl1<0.5 & pdl1reported==1, by(hazard_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 < 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "io_only_test_set_pdl1_under50_xgb_prog_io.png", replace

capture stset progression_days, failure(progression_outcome)
stci, by(hazard_prediction) rmean 
stcox hazard_prediction
sts graph, by(hazard_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Prediction")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "IO Monotherapy Unresponsive" 2 "IO Monotherapy Responsive")) 
graph export "test_set_predictive_io.png", replace

capture stset progression_days if pdl1>=0.5, failure(progression_outcome)
stci if pdl1>=0.5, by(hazard_prediction) rmean 
stcox progressed_prediction if pdl1>=0.5
sts graph if pdl1>=0.5, by(hazard_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Prediction, PD-L1 >= 50%")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "IO Monotherapy Unresponsive" 2 "IO Monotherapy Responsive")) 
graph export "test_set_predictive_over50_io.png", replace

gen new_var = 0
replace new_var =1 if progression_outcome==0
rocreg new_var pdl1 if pdl1reported ==1 
gen new_pred = 1-prog_pred
rocreg new_var new_pred  if pdl1reported ==1 

// Mortality
global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/clf_xgb/"
cd "${path}"
set scheme cleanplots
import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/clf_xgb/test_set_365_10000_1.csv", clear 

gen threshold = 0.132
gen time_limit = 365


gen mort_pred = mort_preds
gen pdl1_over_threshold = (pdl1 >=0.5) 

replace diagnosisyear = 2024 - diagnosisyear
gen mortality_prediction = (mort_pred>=threshold)
gen hazard_prediction  = 0
replace hazard_prediction = 1 if mortality_prediction==0
replace mortality_days = time_limit if mortality_days == 0 | mortality_days >time_limit



capture stset mortality_days, failure(mortality_outcome)
stci, by(hazard_prediction) rmean
stcox hazard_prediction
sts graph, by(hazard_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "test_set_ml_mutations_xgb_mort_io.png", replace


drop if alk==1
drop if egfr==1
drop if ros1==1


capture stset mortality_days if pdl1reported==1, failure(mortality_outcome)
stci if  pdl1reported==1, by(pdl1_over_threshold) rmean
stcox pdl1_over_threshold if  pdl1reported==1
sts graph if  pdl1reported==1, by(pdl1_over_threshold) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by PD-L1 Status; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "PD-L1 <50%" 2 "PD-L1 >=50%"))
graph export "io_test_set_pdl1_mort_io.png", replace

capture stset mortality_days, failure(mortality_outcome)
stci, by(hazard_prediction) rmean
stcox hazard_prediction
sts graph, by(hazard_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "test_set_ml_no_mutations_xgb_mort_io.png", replace

capture stset mortality_days if  pdl1reported==1, failure(mortality_outcome)
stci if  pdl1reported==1, by(hazard_prediction) rmean 
stcox hazard_prediction if  pdl1reported==1
sts graph, by(hazard_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "io_only_test_set_xgb_mort_pdl1reported_io.png", replace


capture stset mortality_days if  pdl1>=0.5, failure(mortality_outcome)
stci if  pdl1>=0.5, by(hazard_prediction) rmean 
stcox hazard_prediction if  pdl1>=0.5
sts graph if  pdl1>=0.5, by(hazard_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "io_only_test_set_pdl1_over50_xgb_mort_io.png", replace

capture stset mortality_days if  pdl1<0.5 & pdl1reported==1, failure(mortality_outcome)
stci if  pdl1<0.5 & pdl1reported==1, by(hazard_prediction) rmean 
stcox hazard_prediction if  pdl1<0.5 & pdl1reported==1
sts graph if  pdl1<0.5 & pdl1reported==1, by(hazard_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 < 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
graph export "io_only_test_set_pdl1_under50_xgb_mort_io.png", replace

capture stset mortality_days, failure(mortality_outcome)
stci, by(hazard_prediction) rmean 
stcox hazard_prediction
sts graph, by(hazard_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Prediction")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "IO Monotherapy Unresponsive" 2 "IO Monotherapy Responsive")) 
graph export "test_set_predictive_io_mort.png", replace

capture stset mortality_days if pdl1>=0.5, failure(mortality_outcome)
stci if pdl1>=0.5, by(hazard_prediction) rmean 
stcox hazard_prediction if pdl1>=0.5
sts graph if pdl1>=0.5, by(hazard_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Prediction, PD-L1 >= 50%")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "IO Monotherapy Unresponsive" 2 "IO Monotherapy Responsive")) 
graph export "test_set_predictive_over50_io_mort.png", replace

gen new_var = 0
replace new_var =1 if mortality_outcome==0
rocreg new_var pdl1 if pdl1reported ==1 
gen new_pred = 1-mort_pred
rocreg new_var new_pred  if pdl1reported ==1 
