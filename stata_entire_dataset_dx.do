/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ib1.race ib1.ecog ib1.histology ib1.stage pdl1 ib1.practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan pdl1_given"
 
 global path "/Users/vahluw/Documents/NSCLC_PDL1/"
 cd "${path}"
 set scheme cleanplots


 import delimited /Users/vahluw/Documents/NSCLC_PDL1/all_data.csv, clear 
 sum age_at_diagnosis
 tab gender
 tab race
 tab stage
 tab histology
 tab smoking_status
 tab pdl1
 tab practice_type
 tab progression_12

 histogram pdl1 if pdl1_given==1, percent bin(10) xtitle("PD-L1 Intensity") color(ebblue)
 graph save "PD_L1_distribution_whole_dataset_include_0", replace


  logistic progression_12  ${indiv_covar}  i.egfr i.braf i.ros1 i.alk i.state pdl1_given i.io_mono i.combo_therapy i.chemo i.other_therapy i.io_mono_used


  logistic progression_12  ${indiv_covar}  i.egfr i.braf i.ros1 i.alk i.state pdl1_given i.io_mono i.combo_therapy i.chemo i.other_therapy i.io_mono_used  if in_test_set==0
 predict logit_progression if in_test_set==1
 rocreg progression_12 logit_progression if in_test_set==1, bseed(20547)

  
   **************************** 
 /* Propensity score matching with all DAG factors --> but this time, kaplan-meier curves show days from treatment initiation until progression*/
 ****************************

 global indiv_covar "ib1.race ib1.ecog ib1.histology ib1.stage pdl1 ib1.practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan pdl1_given"
 
 global path "/Users/vahluw/Documents/NSCLC_PDL1/"
 cd "${path}"
 set scheme cleanplots
 import delimited /Users/vahluw/Documents/NSCLC_PDL1/all_data.csv, clear 
 replace progression_12 = 0 if progression_days == 0
 replace progression_days = 365 if progression_days == 0 | progression_days>365
 
 gen therapy_type = 0 
 replace therapy_type = 1 if io_mono == 1
 replace therapy_type = 2 if combo_therapy == 1
 replace therapy_type = 3 if chemo == 1
 /*
 replace therapy_type = 4 if egfr_drug == 1
 replace therapy_type = 5 if alk_drug == 1
  replace therapy_type = 6 if ros1_drug == 1
 replace therapy_type = 7 if braf_drug == 1
   replace therapy_type = 8 if ras_drug == 1
 replace therapy_type = 9 if secondary_chemo_drug == 1
 */
 
  stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy (Without Matching), Entire Dataset") xtitle ("Survival Time (Days)") ytitle("Proportion at Risk") legend(order(1 "Other Therapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Platinum Chemotherapy")) //5 "EGFR Drug" 6 "ALK Drug" 7 "ROS1 Drug" 8 "BRAF Drug" 9 "RAS Drug" 10 "Secondary Chemotherapy Drug")) 
 graph save "prog_survival_with_mutations_pre_match", replace
 sts test therapy_type, logrank
 
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 //drop if therapy_type>=5 & therapy_type <=9
 
   stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy (Without Matching), Entire Dataset") xtitle ("Survival Time (Days)") ytitle("Proportion at Risk") legend(order(1 "Other Therapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Platinum Chemotherapy"  5 "Secondary Chemotherapy Drug")) 
 graph save "prog_survival_without_mutations_pre_match", replace
 sts test therapy_type, logrank
 
 keep if therapy_type == 1 | therapy_type == 3
 replace therapy_type = 0 if therapy_type==3
 drop if stage ==1 | stage ==18
 

 stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy (Without Matching), Entire Dataset") xtitle ("Survival Time (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph save "prog_survival_chemo_vs_first_line_vs_combo_entire_dataset_zero_tx", replace
 sts test therapy_type, logrank
 
 logit therapy_type ${indiv_covar}
 predict yhat
 
graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching, Entire Dataset") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))
graph save "propensity_score_pre_matching_progression", replace
graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "Chemotherapy") ///
		label(2 "IO Monotherapy"))
graph save "propensity_score_pre_matching_progression_histogram_whole_dataset_tx", replace

pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar}, outcome(progression_12) caliper(0.2) n(1) noreplacement 
count 

sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar}, treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Mono"))
graph save "propensity_score_post_matching_progression_entire_dataset_tx", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "Chemotherapy") ///
		label(2 "IO Monotherapy")) 
graph save "propensity_score_post_matching_progression_histogram_whole_dataset_tx", replace



 stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time From Therapy Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph save "prog_survival_chemo_vs_first_line_entire_dataset_tx", replace
  sts test therapy_type, logrank
  
   stset progression_days if pdl1>=1.5, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset (PD-L1>=50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time From Therapy Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph save "prog_survival_chemo_vs_first_line_entire_dataset_high_pdl1_tx", replace
  sts test therapy_type, logrank
  
 stset progression_days if pdl1<1.5 & pdl1>= 1.0, failure(progression_12)
 stci if pdl1<1.5 & pdl1>= 1.0, by(therapy_type) rmean
 stcox therapy_type if pdl1<1.5 & pdl1>= 1.0
 sts graph if pdl1<1.5 & pdl1>= 1.0, by(therapy_type) title("Progression-Free Survival for Entire Dataset (PD-L1<50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time From Therapy Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph save "prog_survival_chemo_vs_first_line_entire_dataset_low_pdl1_tx", replace
  sts test therapy_type, logrank
  
  