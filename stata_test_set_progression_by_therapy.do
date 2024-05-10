/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */

 global indiv_covar "i.ecog i.histology pdl1 ethnicity i.practice_type diag_year age_at_diagnosis i.race i.gender i.smoking_status days_from_dx_to_tx pt_assistance other_gov_insurance medicare self_pay medicaid commercial_health_plan other_no_insurance i.stage albumin abx steroid creatinine alt ast bilirubin "
 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots

 import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set_365_100.csv", clear 
 rocreg progression_outcome prog_gb_preds
 
  rocreg mortality_outcome mort_gb_preds
 
  import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set_365_100.csv", clear 
 gen threshold = 0.5
 gen time_limit = 365

 
 gen prog_pred = prog_gb_preds
 gen pdl1_over_threshold = (pdl1 >=0.5) 

 replace diag_year = 2024 - diag_year
 gen progressed_prediction = (prog_pred>=threshold)
 
 sum progression_days
 replace progression_days = time_limit if progression_days == 0 | progression_days >time_limit

 
   stset progression_days, failure(progression_outcome)
 stci, by(progressed_prediction) rmean
 stcox progressed_prediction
 sts graph, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_chemo_test_set_ml_mutations_hgb.png", replace

 
 gen therapy_type = -1
 replace therapy_type = 0 if first_line_chemo == 1
 replace therapy_type = 1 if io_mono_used > 0
 replace therapy_type = 2 if combo_therapy == 1
 replace therapy_type = 3 if secondary_chemo_drug == 1
 replace therapy_type = 4 if alk_drug == 1
 replace therapy_type = 5 if egfr_drug == 1
 replace therapy_type = 6 if braf_drug == 1
 replace therapy_type = 7 if ros1_drug == 1
 replace therapy_type = 8  if other_first_line_therapy == 1 | ras_drug == 1
 replace therapy_type = 9 if trk_inhibitor == 1
 replace therapy_type = 10 if met_drug== 1
  replace therapy_type = 11 if carboplatin_only == 1
  replace therapy_type = 12 if cisplatin_only == 1
 

 drop if therapy_type >=2
  
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 
 logit progression_outcome i.therapy_type#c.prog_pred $indiv_covar alk egfr braf  ros1 kras  i.state bev_used three_plus_chemo_drugs  kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year bev_used clinical_study_drug if pdl1_given==1
  logit progression_outcome i.therapy_type#i.progressed_prediction  $indiv_covar alk egfr braf  ros1 kras  i.state bev_used three_plus_chemo_drugs  kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year bev_used clinical_study_drug if pdl1_given==1
  
   logit progression_outcome i.therapy_type#c.pdl1 $indiv_covar alk egfr braf  ros1 kras  i.state bev_used three_plus_chemo_drugs  kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year bev_used clinical_study_drug if pdl1_given==1
  logit progression_outcome i.therapy_type#i.pdl1_over_threshold  $indiv_covar alk egfr braf  ros1 kras  i.state bev_used three_plus_chemo_drugs  kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year bev_used clinical_study_drug if pdl1_given==1

   logit progression_outcome i.therapy_type#c.pdl1 if pdl1_given==1
  logit progression_outcome i.therapy_type#i.pdl1_over_threshold if pdl1_given==1
 
    stset progression_days, failure(progression_outcome)
 stci, by(progressed_prediction) rmean
 stcox progressed_prediction
 sts graph, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_chemo_test_set_ml_no_mutations_hgb.png", replace

 stset progression_days if therapy_type==1, failure(progression_outcome)
 stci if therapy_type==1, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==1
 sts graph if therapy_type==1 , by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_hgb.png", replace
 
  stset progression_days if therapy_type==1 & pdl1>=0.5, failure(progression_outcome)
 stci if therapy_type==1 & pdl1>=0.5, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==1 & pdl1>=0.5
 sts graph if therapy_type==1 & pdl1>=0.5, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_pdl1_over50_hgb.png", replace

   stset progression_days if io_mono_used==2 & pdl1>=0.5, failure(progression_outcome)
 stci if io_mono_used==2 & pdl1>=0.5, by(progressed_prediction) rmean 
 stcox progressed_prediction if io_mono_used==2& pdl1>=0.5
 sts graph if io_mono_used==2 & pdl1>=0.5, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on Pembrolizumab Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "pembro_only_test_set_pdl1_over50_hgb_prog.png", replace
 
 
  stset progression_days if therapy_type==0, failure(progression_outcome)
 stci if therapy_type==0, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==0
 sts graph if therapy_type==0 , by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on Chemotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "chemo_only_test_set_hgb.png", replace
 
 global indiv_covar "ib1.race ib1.ecog ib1.histology ib1.stage pdl1 ib1.practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan"


 gen over_pdl1_threshold = 0
 replace over_pdl1_threshold  = 1 if pdl1>=0.5
 
 
 count if progression_outcome==1 & over_pdl1_threshold ==1 & therapy_type == 1 // Over 50, IO, progressed (incorrect)
 count if progression_outcome==0 & over_pdl1_threshold ==1 & therapy_type == 1 // Over 50, IO, not progressed (correct)
 count if progression_outcome==1 & over_pdl1_threshold ==0 & therapy_type == 1 // Under 50, IO, progressed (correct)
 count if progression_outcome==0 & over_pdl1_threshold ==0 & therapy_type == 1 // Under 50, IO, not progressed (incorrect) 
 count if progression_outcome==1 & over_pdl1_threshold ==1 & therapy_type == 0 // Over 50, chemo, progressed (correct)
 count if progression_outcome==0 & over_pdl1_threshold ==1 & therapy_type == 0 // Over 50, chemo, not progressed (incorrect)
 count if progression_outcome==1 & over_pdl1_threshold ==0 & therapy_type == 0 // Under 50, chemo, progressed (incorrect)
 count if progression_outcome==0 & over_pdl1_threshold ==0 & therapy_type == 0 // Under 50, chemo, not progressed (correct) 
 
 




// Mortality
 import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set_365_100.csv", clear 
  gen threshold = 0.313
 gen prog_pred = mort_gb_preds
 gen pdl1_over_threshold = (pdl1 >=0.5) 
  gen endpoint_prediction = (prog_pred>=threshold)
  gen endpoint = mortality_outcome
 gen time_limit = 365
gen endpoint_days = mortality_days
replace endpoint_days = censor_days if censor_days < mortality_days
replace endpoint_days = censor_days if mortality_days==0 & censor_days < time_limit
 replace diag_year = 2024 - diag_year

 
 sum endpoint_days
 replace endpoint_days = time_limit if endpoint_days == 0 | endpoint_days >time_limit

 
   stset endpoint_days, failure(endpoint)
 stci, by(endpoint_prediction) rmean
 stcox endpoint_prediction
 sts graph, by(endpoint_prediction) title("Overall Survival for Test-Set Patients") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_chemo_test_set_ml_mutations_hgb_overall.png", replace

 gen therapy_type = -1
 replace therapy_type = 0 if first_line_chemo == 1
 replace therapy_type = 1 if io_mono_used > 0
 replace therapy_type = 2 if combo_therapy == 1
 replace therapy_type = 3 if secondary_chemo_drug == 1
 replace therapy_type = 4 if alk_drug == 1
 replace therapy_type = 5 if egfr_drug == 1
 replace therapy_type = 6 if braf_drug == 1
 replace therapy_type = 7 if ros1_drug == 1
 replace therapy_type = 8  if other_first_line_therapy == 1 | ras_drug == 1
 replace therapy_type = 9 if trk_inhibitor == 1
 replace therapy_type = 10 if met_drug== 1
  replace therapy_type = 11 if carboplatin_only == 1
  replace therapy_type = 12 if cisplatin_only == 1

 drop if therapy_type >=2
  
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
  logit mortality_outcome i.therapy_type#c.prog_pred
  logit mortality_outcome i.therapy_type#i.endpoint_prediction
  
   logit mortality_outcome i.therapy_type#c.pdl1 if pdl1_given==1
  logit mortality_outcome i.therapy_type#i.pdl1_over_threshold if pdl1_given==1
 
    stset endpoint_days, failure(endpoint)
 stci, by(endpoint_prediction) rmean
 stcox endpoint_prediction
 sts graph, by(endpoint_prediction) title("Overall Survival for Test-Set Patients") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_chemo_test_set_ml_no_mutations_hgb_overall.png", replace

 stset endpoint_days if therapy_type==1, failure(endpoint)
 stci if therapy_type==1, by(endpoint_prediction) rmean 
 stcox endpoint_prediction if therapy_type==1
 sts graph if therapy_type==1 , by(endpoint_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_hgb_overall.png", replace
 
  stset endpoint_days if therapy_type==1 & pdl1>=0.5, failure(endpoint)
 stci if therapy_type==1 & pdl1>=0.5, by(endpoint_prediction) rmean 
 stcox endpoint_prediction if therapy_type==1 & pdl1>=0.5
 sts graph if therapy_type==1 & pdl1>=0.5, by(endpoint_prediction) title("Overall Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_pdl1_over50_hgb_overall.png", replace

 
   stset endpoint_days if io_mono_used==2 & pdl1>=0.5, failure(endpoint)
 stci if io_mono_used==2 & pdl1>=0.5, by(endpoint_prediction) rmean 
 stcox endpoint_prediction if io_mono_used==2& pdl1>=0.5
 sts graph if io_mono_used==2 & pdl1>=0.5, by(endpoint_prediction) title("Progression-Free Survival for Test-Set Patients on Pembrolizumab Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "pembro_only_test_set_pdl1_over50_mort.png", replace
 
  stset endpoint_days if therapy_type==0, failure(endpoint)
 stci if therapy_type==0, by(endpoint_prediction) rmean 
 stcox endpoint_prediction if therapy_type==0
 sts graph if therapy_type==0 , by(endpoint_prediction) title("Overall Survival for Test-Set Patients on Chemotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Treatment Initiation (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "chemo_only_test_set_hgb_overall.png", replace
 
 global indiv_covar "ib1.race ib1.ecog ib1.histology ib1.stage pdl1 ib1.practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan"


 gen over_pdl1_threshold = 0
 replace over_pdl1_threshold  = 1 if pdl1>=0.5
 
 
 count if progression_outcome==1 & over_pdl1_threshold ==1 & therapy_type == 1 // Over 50, IO, progressed (incorrect)
 count if progression_outcome==0 & over_pdl1_threshold ==1 & therapy_type == 1 // Over 50, IO, not progressed (correct)
 count if progression_outcome==1 & over_pdl1_threshold ==0 & therapy_type == 1 // Under 50, IO, progressed (correct)
 count if progression_outcome==0 & over_pdl1_threshold ==0 & therapy_type == 1 // Under 50, IO, not progressed (incorrect) 
 count if progression_outcome==1 & over_pdl1_threshold ==1 & therapy_type == 0 // Over 50, chemo, progressed (correct)
 count if progression_outcome==0 & over_pdl1_threshold ==1 & therapy_type == 0 // Over 50, chemo, not progressed (incorrect)
 count if progression_outcome==1 & over_pdl1_threshold ==0 & therapy_type == 0 // Under 50, chemo, progressed (incorrect)
 count if progression_outcome==0 & over_pdl1_threshold ==0 & therapy_type == 0 // Under 50, chemo, not progressed (correct) 
 
 
 /* Try to see what's screwing up the model */


 
  import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set_365_100.csv", clear 

 gen threshold = 0.5
 gen prog_pred = prog_gb_preds
 gen pdl1_over_threshold = (pdl1 >=0.5) 
 gen progressed_prediction = (prog_pred>=threshold)
 gen incorrect = 0
 replace incorrect = 1 if (progressed_prediction != progression_outcome)
 
 gen therapy_type = -1
 replace therapy_type = 0 if first_line_chemo == 1
 replace therapy_type = 1 if io_mono_used > 0
 replace therapy_type = 2 if combo_therapy == 1
 replace therapy_type = 3 if secondary_chemo_drug == 1
 replace therapy_type = 4 if alk_drug == 1
 replace therapy_type = 5 if egfr_drug == 1
 replace therapy_type = 6 if braf_drug == 1
 replace therapy_type = 7 if ros1_drug == 1
 replace therapy_type = 8  if other_first_line_therapy == 1 | ras_drug == 1
 replace therapy_type = 9 if trk_inhibitor == 1
 replace therapy_type = 10 if met_drug== 1
  replace therapy_type = 11 if carboplatin_only == 1
  replace therapy_type = 12 if cisplatin_only == 1

 
 logit incorrect i.therapy_type ${indiv_covar} kras  i.state  therapy_year kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets  clinical_study_drug creatinine bilirubin ast alt albumin steroid abx
 
