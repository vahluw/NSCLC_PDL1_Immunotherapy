 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots

 import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set_365_10000.csv", clear 
 rocreg progression_outcome prog_preds
 
  rocreg mortality_outcome mort_preds
 
  import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set_365_10000.csv", clear 
  
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
 
  
 
 gen threshold = 0.704
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

 
   stset progression_days, failure(progression_outcome)
 stci, by(progressed_prediction) rmean
 stcox progressed_prediction
 stcox hazard_prediction
 sts graph, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_chemo_test_set_ml_mutations_xgb_prog.png", replace

 
 gen therapy_type = 1
 replace therapy_type = 0 if firstlinechemotherapy == 1
 replace therapy_type = 2 if firstlinecombinationtherapy == 1
 replace therapy_type = 3 if nonfirstlinechemotherapy == 1
 replace therapy_type = 4 if antialkdrug == 1
 replace therapy_type = 5 if antiegfrdrug == 1
 replace therapy_type = 6 if antibrafdrug == 1
 replace therapy_type = 7 if antiros1drug == 1
 replace therapy_type = 8  if otherfirstlinetherapy == 1 | antirasdrug == 1
 replace therapy_type = 9 if trkinhibitor == 1
 replace therapy_type = 10 if metinhibitor== 1
  replace therapy_type = 11 if carboplatinmonotherapy == 1
  replace therapy_type = 12 if cisplatinmonotherapy == 1
 


  
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 



 stset progression_days if therapy_type==1 & pdl1>=0.01, failure(progression_outcome)
 stci if therapy_type==1 & pdl1>=0.01, by(pdl1_over_threshold) rmean
 stcox pdl1_over_threshold if therapy_type==1 & pdl1>=0.01
 sts graph if therapy_type==1 & pdl1>=0.01, by(pdl1_over_threshold) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by PD-L1 Status; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "1% <= PD-L1 <50%" 2 "PD-L1 >=50%"))
 graph export "io_test_set_pdl1_prog.png", replace
 
  stset progression_days if therapy_type==0 & pdl1>=0.01, failure(progression_outcome)
 stci if therapy_type==0 & pdl1>=0.01, by(pdl1_over_threshold) rmean
 stcox pdl1_over_threshold if therapy_type==0 & pdl1>=0.01
 sts graph if therapy_type==0 & pdl1>=0.01, by(pdl1_over_threshold) title("Progression-Free Survival for Test-Set Patients on Chemotherapy") subtitle("by PD-L1 Status; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "1% <= PD-L1 <50%" 2 "PD-L1 >=50%"))
 graph export "chemo_test_set_pdl1_prog.png", replace
 
   stset progression_days if therapy_type==2 & pdl1>=0.01, failure(progression_outcome)
 stci if therapy_type==2 & pdl1>=0.01, by(pdl1_over_threshold) rmean
 stcox pdl1_over_threshold if therapy_type==2 & pdl1>=0.01
 sts graph if therapy_type==2 & pdl1>=0.01, by(pdl1_over_threshold) title("Progression-Free Survival for Test-Set Patients on Combination Therapy") subtitle("by PD-L1 Status; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "1% <= PD-L1 <50%" 2 "PD-L1 >=50%"))
 graph export "combo_test_set_pdl1_prog.png", replace
 
 
 

 stset mortality_days if therapy_type==1 & pdl1>=0.01, failure(mortality_outcome)
 stci if therapy_type==1 & pdl1>=0.01, by(pdl1_over_threshold) rmean
 stcox pdl1_over_threshold if therapy_type==1 & pdl1>=0.01
 sts graph if therapy_type==1 & pdl1>=0.01, by(pdl1_over_threshold) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by PD-L1 Status; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "1% <= PD-L1 <50%" 2 "PD-L1 >=50%"))
 graph export "io_test_set_pdl1_overall.png", replace
 
  stset mortality_days if therapy_type==0 & pdl1>=0.01, failure(mortality_outcome)
 stci if therapy_type==0 & pdl1>=0.01, by(pdl1_over_threshold) rmean
 stcox pdl1_over_threshold if therapy_type==0 & pdl1>=0.01
 sts graph if therapy_type==0 & pdl1>=0.01, by(pdl1_over_threshold) title("Overall Survival for Test-Set Patients on Chemotherapy") subtitle("by PD-L1 Status; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "1% <= PD-L1 <50%" 2 "PD-L1 >=50%"))
 graph export "chemo_test_set_pdl1_overall.png", replace
 
   stset mortality_days if therapy_type==2 & pdl1>=0.01, failure(mortality_outcome)
 stci if therapy_type==2 & pdl1>=0.01, by(pdl1_over_threshold) rmean
 stcox pdl1_over_threshold if therapy_type==2 & pdl1>=0.01
 sts graph if therapy_type==2 & pdl1>=0.01, by(pdl1_over_threshold) title("Overall Survival for Test-Set Patients on Combination Therapy") subtitle("by PD-L1 Status; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "1% <= PD-L1 <50%" 2 "PD-L1 >=50%"))
 graph export "combo_test_set_pdl1_overall.png", replace
 
 
 
 
drop if therapy_type >=2
  
logit progression_outcome i.therapy_type pdl1 i.therapy_type#c.pdl1 if pdl1reported==1
logit progression_outcome i.therapy_type i.pdl1_over_threshold i.therapy_type#i.pdl1_over_threshold if pdl1reported==1
logit progression_outcome i.therapy_type progressed_prediction i.therapy_type#i.progressed_prediction if pdl1reported==1
 
 

 stset progression_days, failure(progression_outcome)
 stci, by(progressed_prediction) rmean

 stcox progressed_prediction
 stcox hazard_prediction
 sts graph, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients") subtitle("by ML-Derived Risk; Driver-Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_chemo_test_set_ml_no_mutations_xgb_prog.png", replace


 stset progression_days if therapy_type==1, failure(progression_outcome)
 stci if therapy_type==1, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==1
  stcox hazard_prediction if therapy_type==1
 sts graph if therapy_type==1 , by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_xgb_prog.png", replace
 
   teffects ipw (progression_outcome) (progressed_prediction) if therapy_type==1
 
  stset progression_days if therapy_type==1 & pdl1>=0.5, failure(progression_outcome)
 stci if therapy_type==1 & pdl1>=0.5, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==1 & pdl1>=0.5
  stcox hazard_prediction if therapy_type==1 & pdl1>=0.5
 sts graph if therapy_type==1 & pdl1>=0.5, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_pdl1_over50_xgb_prog.png", replace
 
   stset progression_days if therapy_type==0, failure(progression_outcome)
 stci if therapy_type==0, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==0
  stcox hazard_prediction if therapy_type==0
 sts graph if therapy_type==0, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on Chemotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "chemo_only_test_set_pdl1_xgb_prog.png", replace
 
    stset progression_days if firstlinepembrolizumabmonotherap==1, failure(progression_outcome)
 stci if firstlinepembrolizumabmonotherap==1, by(progressed_prediction) rmean 
 stcox progressed_prediction if firstlinepembrolizumabmonotherap==1
  stcox hazard_prediction if firstlinepembrolizumabmonotherap==1
 sts graph if firstlinepembrolizumabmonotherap==1, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on Pembrolizumab Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "pembro_only_test_set_pdl1_xgb_prog.png", replace
 
   
 
 gen stratum = -1
 
 replace stratum = 0 if (therapy_type==0 & progressed_prediction==1)
 replace stratum = 1 if (therapy_type==0 & progressed_prediction==0)
 replace stratum = 2 if (therapy_type==1 & progressed_prediction==1)
 replace stratum = 3 if (therapy_type==1 & progressed_prediction==0)

 
  stset progression_days if therapy_type<=1, failure(progression_outcome)
 stci if therapy_type<=1, by(stratum) rmean 
 stcox progressed_prediction if therapy_type<=1
 sts graph if therapy_type<=1, by(stratum) title("Progression-Free Survival for Test-Set Patients") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Chemo Unresponsive" 2 "Chemo Responsive" 3 "IO Unreponsive" 4 "IO Responsive"))
 graph export "io_chemo_test_set_predictive.png", replace


   stset progression_days if firstlinepembrolizumabmonotherap==1 & pdl1>=0.5, failure(progression_outcome)
 stci if firstlinepembrolizumabmonotherap==1 & pdl1>=0.5, by(progressed_prediction) rmean 
 stcox progressed_prediction if firstlinepembrolizumabmonotherap==1 & pdl1>=0.5
 sts graph if firstlinepembrolizumabmonotherap==1 & pdl1>=0.5, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on Pembrolizumab Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "pembro_only_test_set_pdl1_over50_xgb_prog.png", replace
 


// Mortality
 import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set_365_10000.csv", clear 
  gen threshold = 0.468
 gen prog_pred = mort_preds
 gen pdl1_over_threshold = (pdl1 >=0.5) 
  gen endpoint_prediction = (prog_pred>=threshold)
  gen endpoint = mortality_outcome
 gen time_limit = 365
gen endpoint_days = mortality_days
replace endpoint_days = censor_days if censor_days < mortality_days
replace endpoint_days = censor_days if mortality_days==0 & censor_days < time_limit
 replace diagnosisyear = 2024 - diagnosisyear

 
 sum endpoint_days
 replace endpoint_days = time_limit if endpoint_days == 0 | endpoint_days >time_limit

  gen hazard_prediction  = 0
replace hazard_prediction = 1 if endpoint_prediction==0
   stset endpoint_days, failure(endpoint)
 stci, by(endpoint_prediction) rmean
 stcox endpoint_prediction
 stcox hazard_prediction
 sts graph, by(endpoint_prediction) title("Overall Survival for Test-Set Patients") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_chemo_test_set_ml_mutations_xgb_overall.png", replace

 gen therapy_type = 1
 replace therapy_type = 0 if firstlinechemotherapy == 1
 replace therapy_type = 2 if firstlinecombinationtherapy == 1
 replace therapy_type = 3 if nonfirstlinechemotherapy == 1
 replace therapy_type = 4 if antialkdrug == 1
 replace therapy_type = 5 if antiegfrdrug == 1
 replace therapy_type = 6 if antibrafdrug == 1
 replace therapy_type = 7 if antiros1drug == 1
 replace therapy_type = 8  if otherfirstlinetherapy == 1 | antirasdrug == 1
 replace therapy_type = 9 if trkinhibitor == 1
 replace therapy_type = 10 if metinhibitor== 1
  replace therapy_type = 11 if carboplatinmonotherapy == 1
  replace therapy_type = 12 if cisplatinmonotherapy == 1

 drop if therapy_type >=2
  
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
  logit mortality_outcome i.therapy_type#c.prog_pred
  logit mortality_outcome i.therapy_type#i.endpoint_prediction
  
   logit mortality_outcome i.therapy_type#c.pdl1 if pdl1reported==1
  logit mortality_outcome i.therapy_type#i.pdl1_over_threshold if pdl1reported==1
 
    stset endpoint_days, failure(endpoint)
 stci, by(endpoint_prediction) rmean
 stcox endpoint_prediction
  stcox hazard_prediction
 sts graph, by(endpoint_prediction) title("Overall Survival for Test-Set Patients") subtitle("by ML-Derived Risk; Driver Mutation Negative")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_chemo_test_set_ml_no_mutations_xgb_overall.png", replace

 stset endpoint_days if therapy_type==1, failure(endpoint)
 stci if therapy_type==1, by(endpoint_prediction) rmean 
 stcox endpoint_prediction if therapy_type==1
   stcox hazard_prediction

 sts graph if therapy_type==1 , by(endpoint_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_xgb_overall.png", replace
 
  stset endpoint_days if therapy_type==1 & pdl1>=0.5, failure(endpoint)
 stci if therapy_type==1 & pdl1>=0.5, by(endpoint_prediction) rmean 
 stcox endpoint_prediction if therapy_type==1 & pdl1>=0.5
   stcox hazard_prediction

 sts graph if therapy_type==1 & pdl1>=0.5, by(endpoint_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_pdl1_over50_xgb_overall.png", replace

 
   stset endpoint_days if firstlinepembrolizumabmonotherap==1 & pdl1>=0.5, failure(endpoint)
 stci if firstlinepembrolizumabmonotherap==1& pdl1>=0.5, by(endpoint_prediction) rmean 
 stcox endpoint_prediction if firstlinepembrolizumabmonotherap==1 & pdl1>=0.5
   stcox hazard_prediction

 sts graph if firstlinepembrolizumabmonotherap==1 & pdl1>=0.5, by(endpoint_prediction) title("Overall Survival for Test-Set Patients on Pembrolizumab Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "pembro_only_test_set_pdl1_over50_overall.png", replace
 
  stset endpoint_days if therapy_type==0, failure(endpoint)
 stci if therapy_type==0, by(endpoint_prediction) rmean 
 stcox endpoint_prediction if therapy_type==0
   stcox hazard_prediction

 sts graph if therapy_type==0 , by(endpoint_prediction) title("Overall Survival for Test-Set Patients on Chemotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "chemo_only_test_set_xgb_overall.png", replace



 
  import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set_365_10000.csv", clear 

 gen threshold = 0.66
 gen prog_pred = prog_preds
 gen pdl1_over_threshold = (pdl1 >=0.5) 
 gen progressed_prediction = (prog_pred>=threshold)
 gen incorrect = 0
 replace incorrect = 1 if (progressed_prediction != progression_outcome)
 
 gen therapy_type = 1
 replace therapy_type = 0 if firstlinechemotherapy == 1
 //replace therapy_type = 1 if io_mono_used > 0
 replace therapy_type = 2 if firstlinecombinationtherapy == 1
 replace therapy_type = 3 if nonfirstlinechemotherapy == 1
 replace therapy_type = 4 if antialkdrug == 1
 replace therapy_type = 5 if antiegfrdrug == 1
 replace therapy_type = 6 if antibrafdrug == 1
 replace therapy_type = 7 if antiros1drug == 1
 replace therapy_type = 8  if otherfirstlinetherapy == 1 | antirasdrug == 1
 replace therapy_type = 9 if trkinhibitor == 1
 replace therapy_type = 10 if metinhibitor== 1
  replace therapy_type = 11 if carboplatinmonotherapy == 1
  replace therapy_type = 12 if cisplatinmonotherapy == 1

 
 logit incorrect i.therapy_type  kras  therapyyear  chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecificedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen
 
