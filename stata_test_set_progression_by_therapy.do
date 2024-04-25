/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */


 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots

  
 import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/test_set.csv", clear 
 rocreg progression_12 gb_preds, bseed(0)
 gen threshold = 0.37
 
 gen prog_pred = gb_preds
 gen pdl1_over_threshold = (pdl1 >=0.5) 

 replace diag_year = 2024 - diag_year
 gen progressed_prediction = (prog_pred>=threshold)
 
 sum progression_days
 replace progression_days = 366 if progression_days == 0 | progression_days >365 

 
   stset progression_days, failure(progression_12)
 stci, by(progressed_prediction) rmean
 stcox progressed_prediction
 sts graph, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients") subtitle("by ML-Derived Risk, Pre-Matching")  xtitle ("Survival Time From Diagnosis Date (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "prog_survival_chemo_vs_first_line_test_high_risk_vs_low_risk_test_set.png", replace
 
 /*
 io_mono, io_mono_used, combo_therapy, first_line_chemo, secondary_chemo_drug, other_therapy,
                        alk_drug, egfr_drug, braf_drug, ros1_drug, ras_drug, other_first_line_therapy,
 */
 gen therapy_type = 11
 replace therapy_type = 0 if first_line_chemo == 1
 replace therapy_type = 1 if io_mono == 1
 replace therapy_type = 2 if combo_therapy == 1
 replace therapy_type = 3 if secondary_chemo_drug == 1
 replace therapy_type = 4 if other_therapy == 1
 replace therapy_type = 5 if alk_drug == 1
 replace therapy_type = 6 if egfr_drug == 1
 replace therapy_type = 7 if braf_drug == 1
 replace therapy_type = 8 if ros1_drug == 1
 replace therapy_type = 9 if ras_drug == 1
 replace therapy_type = 10 if other_first_line_therapy == 1


 drop if therapy_type >=2
  
 drop if alk==1
 drop if egfr==1
 drop if ros1==1

 stset progression_days if therapy_type==1, failure(progression_12)
 stci if therapy_type==1, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==1
 sts graph if therapy_type==1 , by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Diagnosis Date  (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set.png", replace
 
  stset progression_days if therapy_type==1 & pdl1>=0.5, failure(progression_12)
 stci if therapy_type==1 & pdl1>=0.5, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==1 & pdl1>=0.5
 sts graph if therapy_type==1 & pdl1>=0.5, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, PD-L1 >= 50% Only")  xtitle ("Survival Time From Diagnosis Date  (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_set_pdl1_over50.png", replace

 
  stset progression_days if therapy_type==0, failure(progression_12)
 stci if therapy_type==0, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==0
 sts graph if therapy_type==0 , by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on Doublet Chemotherapy") subtitle("by ML-Derived Risk")  xtitle ("Survival Time From Diagnosis Date  (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "chemo_only_test_set.png", replace
 
 global indiv_covar "ib1.race ib1.ecog ib1.histology ib1.stage pdl1 ib1.practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan"

 keep if pdl1_given==1
 gen over_pdl1_threshold = 0
 replace over_pdl1_threshold  = 1 if pdl1>=0.5
 
 
 count if progression_12==1 & over_pdl1_threshold ==1 & therapy_type == 1 // Over 50, IO, progressed (incorrect)
 count if progression_12==0 & over_pdl1_threshold ==1 & therapy_type == 1 // Over 50, IO, not progressed (correct)
 count if progression_12==1 & over_pdl1_threshold ==0 & therapy_type == 1 // Under 50, IO, progressed (correct)
 count if progression_12==0 & over_pdl1_threshold ==0 & therapy_type == 1 // Under 50, IO, not progressed (incorrect) 
 count if progression_12==1 & over_pdl1_threshold ==1 & therapy_type == 0 // Over 50, chemo, progressed (correct)
 count if progression_12==0 & over_pdl1_threshold ==1 & therapy_type == 0 // Over 50, chemo, not progressed (incorrect)
 count if progression_12==1 & over_pdl1_threshold ==0 & therapy_type == 0 // Under 50, chemo, progressed (incorrect)
 count if progression_12==0 & over_pdl1_threshold ==0 & therapy_type == 0 // Under 50, chemo, not progressed (correct) 
 
 


 