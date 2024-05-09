   **************************** 
 /* Propensity score matching with all DAG factors --> but this time, kaplan-meier curves show days from treatment initiation until progression*/
 ****************************

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "i.ecog i.histology pdl1 ethnicity i.practice_type diag_year age_at_diagnosis i.race i.gender i.smoking_status days_from_dx_to_tx pt_assistance other_gov_insurance medicare self_pay medicaid commercial_health_plan other_no_insurance i.stage albumin abx steroid creatinine alt ast bilirubin "
 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots
 



 import delimited "all_data_365_100.csv", clear 

  gen time_limit = 365
 gen outcome = "progression"

 gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5
 
 gen censor_time = progression_days 
 replace censor_time = mortality_days if outcome == "mortality"
 gen endpoint = progression_outcome
 replace endpoint = mortality_outcome if outcome == "mortality"
 
 sum age_at_diagnosis
 tab gender
 tab race
 tab stage
 tab histology
 tab smoking_status
 tab pdl1
 tab practice_type
 tab progression_outcome
 gen insured = 0
 replace insured = 1 if pt_assistance == 1 | other_gov_insurance == 1 | medicare == 1 | medicaid == 1 | commercial_health_plan ==1


 histogram pdl1, percent bin(9) xtitle("PD-L1 Intensity") color(ebblue)
 graph export "PD_L1_distribution_whole_dataset_include_0.png", replace
 
 
 drop if ras_drug==1
 gen therapy_type = -1
 replace therapy_type = 0 if first_line_chemo == 1
 replace therapy_type = 1 if io_mono_used > 0
 replace therapy_type = 2 if combo_therapy == 1
 replace therapy_type = 3 if secondary_chemo_drug == 1
 replace therapy_type = 4 if alk_drug == 1
 replace therapy_type = 5 if egfr_drug == 1
 replace therapy_type = 6 if braf_drug == 1
 replace therapy_type = 7 if ros1_drug == 1
 replace therapy_type = 8  if other_first_line_therapy == 1 
 replace therapy_type = 9 if trk_inhibitor == 1
 replace therapy_type = 10 if met_drug== 1
  replace therapy_type = 11 if carboplatin_only == 1
  replace therapy_type = 12 if cisplatin_only == 1
 

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 logit endpoint i.therapy_type ${indiv_covar} alk egfr braf  ros1 kras  i.state bev_used three_plus_chemo_drugs  kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year bev_used clinical_study_drug
 predict logit_pred
 rocreg endpoint logit_pred
 
 
  stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Platinum-Based Chemotherapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Other Chemotherapy" 5 "ALK Drug" 6 "EGFR Drug" 7 "BRAF Drug" 8 "ROS1 Drug" 9  "Other First-Line Therapy" 10 "TRK Inhibitor" 11 "MET Drug" 12 "Carboplatin Only" 13 "Cisplatin Only"))
 graph export "prog_survival_with_mutations_pre_match_prog.png", replace
 sts test therapy_type, logrank
 
 

 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if therapy_type>=2
 drop if bev_used == 1 | three_plus_chemo_drugs == 1
 drop if clinical_study_drug == 1
 
 drop if stage == 1 | stage == 18
 drop if race == 4  | stage == 4
 teffects psmatch (progression_outcome) (therapy_type $indiv_covar pdl1_given)
teffects ipwra (progression_outcome $indiv_covar i.state bev_used three_plus_chemo_drugs  kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year bev_used clinical_study_drug, logit) (therapy_type $indiv_covar pdl1_given i.state bev_used three_plus_chemo_drugs  kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year bev_used clinical_study_drug)
drop if stage ==5
/*
 teffects psmatch (progression_outcome) (therapy_type $indiv_covar) if pdl1>=0.5
teffects ipwra (progression_outcome $indiv_covar, logit) (therapy_type $indiv_covar) if pdl1>=0.5
 teffects psmatch (progression_outcome) (therapy_type $indiv_covar) if pdl1<0.5 & pdl1>0
 drop if stage==14 | ecog ==5
teffects ipwra (progression_outcome $indiv_covar, logit) (therapy_type $indiv_covar) if pdl1<0.5 & pdl1>0
*/
 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation  (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "prog_survival_without_mutations_pre_match_no_combo_prog.png", replace
 sts test therapy_type, logrank


 logit therapy_type ${indiv_covar} braf kras kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis  interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year  //i.therapy_type##over_threshold
 
predict yhat
graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))
graph export "propensity_score_pre_matching_progression.png", replace
graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))
graph export "propensity_pre_match_hist_prog.png", replace

pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} pdl1_given braf kras kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis  interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year, outcome(endpoint) caliper(0.2) //n(1) noreplacement 
count 


sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))
graph export "propensity_post_match.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 
graph export "propensity_post_match_hist_prog.png", replace



 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_prog.png", replace
  sts test therapy_type, logrank
  
   stset censor_time if pdl1>=0.5, failure(endpoint)
 stci if pdl1>=0.5, by(therapy_type) rmean
 stcox therapy_type if pdl1>=0.5
 sts graph if pdl1>=0.5, by(therapy_type) title("Progression-Free Survival for Entire Dataset (PD-L1>=50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_pdl1_over50_prog.png", replace
  sts test therapy_type, logrank
  
     stset censor_time if pdl1==0.0 & pdl1_given==1, failure(endpoint)
 stci if pdl1==0.0 & pdl1_given==1, by(therapy_type) rmean
 stcox therapy_type if pdl1==0.0 & pdl1_given==1
 sts graph if pdl1==0.0 & pdl1_given==1, by(therapy_type) title("Progression-Free Survival for Entire Dataset (Confirmed PD-L1 Negative)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_pdl1_zero_prog.png", replace
  sts test therapy_type, logrank
  
 keep if pdl1>=0.01 & pdl1 < 0.5
 stset censor_time , failure(endpoint)
 stci , by(therapy_type) rmean
 stcox therapy_type 
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset (1%<=PD-L1<50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_pdl1_under50_prog.png", replace
  sts test therapy_type, logrank
  

  
  
     **************************** 
 /* Propensity score matching with all DAG factors --> but this time, kaplan-meier curves show days from treatment initiation until progression*/
 ****************************

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "i.ecog i.histology pdl1 ethnicity i.practice_type diag_year age_at_diagnosis i.race i.gender i.smoking_status days_from_dx_to_tx pt_assistance other_gov_insurance medicare self_pay medicaid commercial_health_plan other_no_insurance i.stage albumin abx steroid creatinine alt ast bilirubin"
 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots


 import delimited "all_data_365_100.csv", clear 

  gen time_limit = 365
 gen outcome = "mortality"

 
 
 gen censor_time = progression_days 
 replace censor_time = mortality_days if outcome == "mortality"
 gen endpoint = progression_outcome
 replace endpoint = mortality_outcome if outcome == "mortality"
 
 sum age_at_diagnosis
 tab gender
 tab race
 tab stage
 tab histology
 tab smoking_status
 tab pdl1
 tab practice_type
 tab progression_outcome
 gen insured = 0
 replace insured = 1 if pt_assistance == 1 | other_gov_insurance == 1 | medicare == 1 | medicaid == 1 | commercial_health_plan ==1


 histogram pdl1, percent bin(9) xtitle("PD-L1 Intensity") color(ebblue)
 graph export "PD_L1_distribution_whole_dataset_include_0.png", replace
 
 
 
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
 

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 logit endpoint i.therapy_type ${indiv_covar} alk egfr braf  ros1 kras  i.state bev_used three_plus_chemo_drugs  kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year bev_used clinical_study_drug
 predict logit_pred
 rocreg endpoint logit_pred
 
 
  stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Platinum-Based Chemotherapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Other Chemotherapy" 5 "ALK Drug" 6 "EGFR Drug" 7 "BRAF Drug" 8 "ROS1 Drug" 9  "Other First-Line Therapy" 10 "TRK Inhibitor" 11 "MET Drug" 12 "Carboplatin Only" 13 "Cisplatin Only"))
 graph export "overall_survival_with_mutations_pre_match_mortality.png", replace
 sts test therapy_type, logrank
 
 

 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if therapy_type>=2
 drop if bev_used == 1 | three_plus_chemo_drugs == 1
 drop if clinical_study_drug == 1
 
 //drop if stage == 1 | stage == 18
 
 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation  (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "overall_survival_without_mutations_pre_match_no_combo_mortality.png", replace
 sts test therapy_type, logrank

 logit therapy_type ${indiv_covar} braf kras kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets therapy_year  therapy_type*pdl1
 
predict yhat
graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))
graph export "propensity_score_pre_matching_progression.png", replace
graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))
graph export "propensity_pre_match_hist_mortality.png", replace

pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} , outcome(endpoint) caliper(0.2) //n(1) noreplacement 
count 

sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Mono"))
graph export "propensity_post_match.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 
graph export "propensity_post_match_hist_mortality.png", replace



 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for Entire Dataset") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "overall_chemo_io_post_match_mortality.png", replace
  sts test therapy_type, logrank
  
   stset censor_time if pdl1>=0.5, failure(endpoint)
 stci if pdl1>=0.5, by(therapy_type) rmean
 stcox therapy_type if pdl1>=0.5
 sts graph if pdl1>=0.5, by(therapy_type) title("Overall Survival for Entire Dataset (PD-L1>=50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "overall_chemo_io_post_match_pdl1_over50_mortality.png", replace
  sts test therapy_type, logrank
  
     stset censor_time if pdl1==0.0 & pdl1_given==1, failure(endpoint)
 stci if pdl1==0.0 & pdl1_given==1, by(therapy_type) rmean
 stcox therapy_type if pdl1==0.0 & pdl1_given==1
 sts graph if pdl1==0.0 & pdl1_given==1, by(therapy_type) title("Overall Survival for Entire Dataset (Confirmed PD-L1 Negative)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "overall_chemo_io_post_match_pdl1_zero_mortality.png", replace
  sts test therapy_type, logrank
  
 keep if pdl1>=0.01 & pdl1 < 0.5
 stset censor_time , failure(endpoint)
 stci , by(therapy_type) rmean
 stcox therapy_type 
 sts graph, by(therapy_type) title("Overall Survival for Entire Dataset (1%<=PD-L1<50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "overal_chemo_io_post_match_pdl1_under50_mortality.png", replace
  sts test therapy_type, logrank
 
 ////// Regression discontinuity ///////////

import delimited "all_data_365_100.csv", clear
gen therapy_type = -1
replace therapy_type = 0 if first_line_chemo == 1
replace therapy_type = 1 if io_mono >0
drop if io_mono_used == 1 | io_mono_used==6
//keep if io_mono_used==2
keep if therapy_type>=0
gen first_line = 0
replace first_line = 1 if therapy_type ==1
gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5
drop if alk==1
drop if egfr==1
drop if ros1==1
drop if three_plus_chemo_drugs == 1
drop if bev_used==1
keep if pdl1> 0
keep if pdl1_given==1

gen black = 0
gen white = 0
gen asian = 0
gen other_race = 0
gen hispanic = 0
replace black =1 if race== 5
replace white = 1 if race == 1
replace asian = 1 if race == 2
replace other_race = 1 if race == 3
replace hispanic = 1 if ethnicity == 1
gen male = 0
replace male =1 if gender == 2
gen female = 0
replace female = 1 if gender == 1
gen never_smoker = 0
replace never_smoker = 1 if smoking_status==1
gen prev_smoker = 0
replace prev_smoker =1 if smoking_status==2


gen insured = 0
replace insured = 1 if pt_assistance == 1 | other_gov_insurance == 1 | medicare == 1 | medicaid == 1 | commercial_health_plan == 1 | self_pay == 1 

rdplot first_line pdl1, c(0.5) 
graph export "polynomial_fit_RD.png", replace
binscatter first_line pdl1, rd(0.5) yti("Proportion Receiving IO Monotherapy Treatment") xti("PD-L1") 
graph export "discontinuity_treat.png", replace 
binscatter days_from_dx_to_tx pdl1, rd(0.5) yti("Days from Dx to Tx") xti("PD-L1") 
graph export "discontinuity_days.png", replace 
binscatter therapy_year pdl1, rd(0.5) yti("Therapy Year") xti("PD-L1") 
graph export "discontinuity_year.png", replace 
binscatter age_at_diagnosis pdl1, rd(0.5) yti("Age at Diagnosis") xti("PD-L1") 
graph export "discontinuity_age.png", replace 
binscatter practice_type pdl1, rd(0.5) yti("Proportion at Community Center") xti("PD-L1") 
graph export "discontinuity_practice_type.png", replace 
binscatter insured pdl1, rd(0.5) yti("Probability Insured") xti("PD-L1") 
graph export "discontinuity_insured.png", replace 
binscatter never_smoker pdl1, rd(0.5) yti("Proportion Never Smoker") xti("PD-L1") 
graph export "discontinuity_never_smoker.png", replace
binscatter prev_smoker pdl1, rd(0.5) yti("Proportion Previous Smoker") xti("PD-L1") 
graph export "discontinuity_prev_smoker.png", replace
binscatter male pdl1, rd(0.5) yti("Proportion Male") xti("PD-L1") 
graph export "discontinuity_male.png", replace
 
kdensity pdl1 , xline(0.5)

//Plotting all, testing only within the optimal bandwidth estimated
rddensity pdl1 , pl c(0.5)


// Note: The "T" is your local average treatment effect. The P>|T| in the large sample is your p-value. Use the confidence interval.
// RD density plot shows that while there is a significant spike from 0.4 to 0.5 there is no reason to think pathologists
// are artificially inflating PD-L1 values so as to increase likelihood of IO monotherapy. 
/* Actually do statistical analysis for RD without nivolumab for IO vs chemo */
rddensity pdl1, c(0.5) vce(jackknife) plot

rdwinselect pdl1 days_from_dx_to_tx therapy_year age_at_diagnosis practice_type insured black white asian other_race hispanic male female never_smoker prev_smoker, c(0.5) seed(0) reps(1000) level(0.05) wmass
rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.405) wr(0.627)  firststage // Unadjusted
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.0) wr(1) wmass firststage  // Unadjusted
rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.405) wr(0.627) wmass firststage // Adjusted
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.4) wr(0.501) wmass firststage  // Adjusted

drop if therapy_year < 2014

rdplot first_line pdl1, c(0.5) ci(95) p(3)
graph export "polynomial_fit_RD_after_2014.png", replace
binscatter first_line pdl1, rd(0.5) yti("Probability of IO Monotherapy Treatment") xti("PD-L1") 
graph export "discontinuity_treat_after_2014.png", replace 
binscatter progression_outcome pdl1, rd(0.5) yti("Proportion Experiencing Disease Progression") xti("PD-L1") 
graph export "discontinuity_outcome_after_2014.png", replace 

graph twoway (hist pdl1) , xline(0.5, lcolor(red))
kdensity pdl1 , xline(0.5)

rddensity pdl1 , pl c(0.5)

rddensity pdl1, c(0.5) vce(jackknife) plot
rdwinselect pdl1 days_from_dx_to_tx therapy_year age_at_diagnosis practice_type insured black white asian other_race hispanic male female never_smoker prev_smoker, c(0.5) seed(0) reps(1000)  wmass level(0.15) 
rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.4) wr(0.501) wmass firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.4) wr(0.501) wmass firststage 
  /////////////////////////////
/*  Instrumental variables */

global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
cd "${path}"
set scheme cleanplots
global indiv_covar "i.ecog i.histology  pdl1 ethnicity i.practice_type diag_year age_at_diagnosis i.race i.gender i.smoking_status days_from_dx_to_tx pt_assistance other_gov_insurance medicare medicaid commercial_health_plan other_no_insurance kras braf pdl1_given albumin abx steroid pembrolizumab_used"
  
import delimited "all_data_365_100.csv", clear 

gen therapy_type = -1
replace therapy_type = 0 if first_line_chemo == 1
replace therapy_type = 1 if io_mono_used > 0
drop if alk==1
drop if egfr==1
drop if ros1==1
drop if bev_used == 1 | three_plus_chemo_drugs == 1
keep if therapy_type >= 0
set emptycells drop

logit progression_outcome i.therapy_type $indiv_covar  i.stage
gen over_threshold =(pdl1>=0.5)
//keep if pdl1_given==1
ivreg2 progression_outcome (therapy_type=i.practiceid) // Unadjusted regression
ivreg2 progression_outcome (therapy_type=i.practiceid ) $indiv_covar  i.state  therapy_year kidney_failure chronic_kidney_disease renal_disease kidney_transplant cirrhosis hepatitis liver_transplant connective_tissue scleroderma lupus rheumatoid_arthritis  interstitial_lung_disease diabetes bone_mets brain_mets cns_mets digestive_mets adrenal_mets unspecified_mets ast alt creatinine bilirubin // Adjusted regression

