/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */

 global path "/Users/vahluw/Documents/NSCLC_PDL1/"
 cd "${path}"
 set scheme cleanplots
  import delimited "all_data.csv", clear 
  
 import delimited "test_set.csv", clear 
 rocreg progression_12 ml_static_preds, bseed(0)
 gen threshold = 0.576601
 
 gen prog_pred = ml_static_preds
 gen pdl1_given = pdl1>=1.0
 gen pdl1_over_threshold = (pdl1 >=1.5) 

 replace diag_year = 2024 - diag_year
 gen progressed_prediction = (prog_pred>=threshold)
 
 sum progression_days
 replace progression_days = 366 if progression_days == 0 | progression_days >365 
 
 drop if alk==2
 drop if egfr==2
 drop if ros1==2
 
   stset progression_days, failure(progression_12)
 stci, by(progressed_prediction) rmean
 stcox progressed_prediction
 sts graph, by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients") subtitle("by ML-Derived Risk, Pre-Matching")  xtitle ("Survival Time From Diagnosis Date (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "prog_survival_chemo_vs_first_line_test_high_risk_vs_low_risk_test_set.png", replace
 
 gen therapy_type = 3
 replace therapy_type = 0 if chemo == 1
 replace therapy_type = 1 if first_line == 1
 replace therapy_type = 2 if combo == 1 

  stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean 
 stcox progressed_prediction
 sts graph , by(therapy_type) title("Progression-Free Survival for Test-Set Patients") subtitle("by ML-Derived Risk, Pre-Matching")  xtitle ("Survival Time From Diagnosis Date  (Days)") ytitle ("Proportion at Risk") legend(order(1 "Doublet Chemotherapy" 2 "IO Monotherapy" 3 "Combination Therapy" 4 "Other Therapy"))
 graph export "all_therapies_test_set_predicted.png", replace
 
 drop if therapy_type == 2 | therapy_type == 3
 

 stset progression_days if therapy_type==1, failure(progression_12)
 stci if therapy_type==1, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==1
 sts graph if therapy_type==1 , by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on IO Monotherapy") subtitle("by ML-Derived Risk, Pre-Matching")  xtitle ("Survival Time From Diagnosis Date  (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "io_only_test_high_risk_vs_low_risk_test_set.png", replace
 
  stset progression_days if therapy_type==0, failure(progression_12)
 stci if therapy_type==0, by(progressed_prediction) rmean 
 stcox progressed_prediction if therapy_type==0
 sts graph if therapy_type==0 , by(progressed_prediction) title("Progression-Free Survival for Test-Set Patients on Doublet Chemotherapy") subtitle("by ML-Derived Risk, Pre-Matching")  xtitle ("Survival Time From Diagnosis Date  (Days)") ytitle ("Proportion at Risk") legend(order(1 "Low-Risk" 2 "High-Risk"))
 graph export "chemo_only_test_high_risk_vs_low_risk_test_set.png", replace
 
 global indiv_covar "ib1.race ib1.ecog ib1.histology ib1.stage pdl1 ib1.practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan"

 keep if pdl1>=1.0
 replace pdl1=pdl1-1.0
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
 
 
 
logit over_pdl1_threshold ${indiv_covar}
predict yhat


graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "Chemo") label(2 "IO Mono"))
graph export "propensity_score_pre_matching_progression.png", replace
graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title(Distribution of Propensity Scores by Therapy Type (Pre-Matching)) ///
		subtitle(Test Set) legend(label(1 "Chemotherapy") ///
		label(2 "IO Monotherapy"))
graph export "propensity_score_pre_matching_progression_histogram_test_set.png", replace

pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar}, outcome(progression_12) caliper(0.2) n(1)  noreplacement 

drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar}, treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Post-Matching") xtitle("Propensity Score") legend(label (1 "Chemo") label(2 "IO Monotherapy"))
graph export "propensity_score_post_matching_progression_test_set.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title(Distribution of Propensity Scores by Therapy Type (Post-Matching)) ///
		subtitle(Test Set) legend(label(1 "Chemotherapy") ///
		label(2 "IO Monotherapy"))
graph export "propensity_score_post_matching_progression_histogram_test_set.png", replace


stset progression_days if prog_pred<threshold & pdl1_given==1, failure(progression_12)
 stci if prog_pred<threshold & pdl1_given==1, by(pdl1_over_threshold) rmean
 stcox pdl1_over_threshold if prog_pred<threshold & pdl1_given==1
 sts graph if prog_pred<threshold & pdl1_given==1, by(pdl1_over_threshold) title("Progression-Free Survival for 'Low-Risk' Test-Set Patients") subtitle("by PDL1 Status, Post-Matching")  xtitle ("Survival Time From Diagnosis Date (Days)") ytitle ("Proportion at Risk") legend(order(1 "PDL1 < 50%" 2 "PDL1 >= 50%"))
 graph export "prog_survival_pdl1_test_set_low_risk.png", replace
 
 stset progression_days if prog_pred>=threshold & pdl1_given==1, failure(progression_12)
 stci if prog_pred>=threshold & pdl1_given==1, by(pdl1_over_threshold) rmean
 stcox pdl1_over_threshold  if prog_pred>=threshold & pdl1_given==1
 sts graph if prog_pred>=threshold & pdl1_given==1, by(pdl1_over_threshold) title("Progression-Free Survival for 'High-Risk' Test-Set Patients") subtitle("by PDL1 Status, Post-Matching")  xtitle ("Survival Time From Diagnosis Date (Days)") ytitle ("Proportion at Risk") legend(order(1 "PDL1 < 50%" 2 "PDL1 >= 50%"))
 graph export "prog_survival_pdl1_test_set_high_risk.png", replace
 

 

 