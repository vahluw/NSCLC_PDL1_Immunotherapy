/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "i.race i.ecog i.histology i.stage pdl1 i.practice_type diag_year age_at_diagnosis pdl1_given"


 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots

  
 import delimited "all_data_365.csv", clear 
 sum age_at_diagnosis
 tab gender
 tab race
 tab stage
 tab histology
 tab smoking_status
 tab pdl1
 tab practice_type
 tab progression_12
  gen insured = 0
 replace insured = 1 if pt_assistance == 1 | other_gov_insurance == 1 | medicare == 1 | medicaid == 1 | commercial_health_plan ==1

 histogram pdl1 if pdl1_given==1, percent bin(9) xtitle("PD-L1 Intensity") color(ebblue)
 graph export "PD_L1_distribution_whole_dataset_include_0.png", replace
  
   **************************** 
 /* Propensity score matching with all DAG factors --> but this time, kaplan-meier curves show days from treatment initiation until progression*/
 ****************************

 global indiv_covar "i.race i.ecog i.histology i.stage pdl1 i.practice_type diag_year age_at_diagnosis  pdl1_given"

 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots
 import delimited all_data_365.csv, clear 
 replace progression_12 = 0 if progression_days == 0
 replace progression_days = 365 if progression_days == 0 | progression_days>365
  gen insured = 0
 replace insured = 1 if pt_assistance == 1 | other_gov_insurance == 1 | medicare == 1 | medicaid == 1 | commercial_health_plan ==1
 gen therapy_type = -1
 replace therapy_type = 0 if first_line_chemo == 1
 replace therapy_type = 1 if io_mono == 1
 replace therapy_type = 2 if combo_therapy == 1
 replace therapy_type = 3 if secondary_chemo_drug == 1
 replace therapy_type = 4 if alk_drug == 1
 replace therapy_type = 5 if egfr_drug == 1
 replace therapy_type = 6 if braf_drug == 1
 replace therapy_type = 7 if ros1_drug == 1
 replace therapy_type = 8 if ras_drug == 1
 replace therapy_type = 9 if other_first_line_therapy == 1

 
  stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time (Days)") ytitle("Proportion at Risk") legend(order(1 "Platinum-Based Chemotherapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Other Chemotherapy" 5 "ALK Drug" 6 "EGFR Drug" 7 "BRAF Drug" 8 "ROS1 Drug" 9 "RAS Drug" 10 "Other First-Line Therapy"))
 graph export "prog_survival_with_mutations_pre_match.png", replace
 sts test therapy_type, logrank
 
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if therapy_type>2

 
   stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemotherapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy")) 
 graph export "prog_survival_without_mutations_pre_match_combo.png", replace
 sts test therapy_type, logrank
 drop if therapy_type==2
 
 
   stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time (Days)") ytitle("Proportion at Risk") legend(order(1 "Other Therapy" 2 "First-Line IO Monotherapy")) 
 graph export "prog_survival_without_mutations_pre_match_no_combo.png", replace
 sts test therapy_type, logrank
 
 
logit therapy_type ${indiv_covar} insured
predict yhat
drop if yhat < 0.00001
drop if stage == 1 | stage == 18
teffects ipwra (progression_12 ${indiv_covar} insured pdl1_given) (therapy_type ${indiv_covar} pdl1_given insured)
teffects ipw (progression_12 ) (therapy_type ${indiv_covar} insured pdl1_given)

drop if race==4 | stage==4

teffects ipwra (progression_12 ${indiv_covar} insured) (therapy_type ${indiv_covar} insured) if pdl1>=0.5
teffects ipw (progression_12 ) (therapy_type ${indiv_covar} insured) if pdl1>=0.5


teffects ipwra (progression_12 ${indiv_covar} insured) (therapy_type ${indiv_covar} insured) if pdl1<0.5 & pdl1_given==1
teffects ipw (progression_12 ) (therapy_type ${indiv_covar} insured) if pdl1<0.5 & pdl1_given==1

gen progression_6_months = 0
replace progression_6_months = 1 if progression_days <=182

 
graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))
graph export "propensity_score_pre_matching_progression.png", replace
graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "Chemotherapy") ///
		label(2 "IO Monotherapy"))
graph export "propensity_pre_match_hist.png", replace

pstest ${indiv_covar} insured pdl1_given, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar}  insured pdl1_given, outcome(progression_12) caliper(0.2) n(1) noreplacement 
count 

sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} pdl1_given insured, treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Mono"))
graph export "propensity_post_match.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "Chemotherapy") ///
		label(2 "IO Monotherapy")) 
graph export "propensity_post_match_hist.png", replace



 stset progression_days, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time From Therapy Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match.png", replace
  sts test therapy_type, logrank
  
   stset progression_days if pdl1>=0.5, failure(progression_12)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset (PD-L1>=50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time From Therapy Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_pdl1_over50.png", replace
  sts test therapy_type, logrank
  
 stset progression_days if pdl1<=0.5 & pdl1_given==1, failure(progression_12)
 stci if pdl1<=0.5 & pdl1_given==1, by(therapy_type) rmean
 stcox therapy_type if pdl1<=0.5 & pdl1_given==1
 sts graph if pdl1<=0.5 & pdl1_given==1, by(therapy_type) title("Progression-Free Survival for Entire Dataset (PD-L1<50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time From Therapy Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_pdl1_under50.png", replace
  sts test therapy_type, logrank
  
  /////////////////////////////
/*  Instrumental variables */
 global indiv_covar "i.race i.ecog i.histology i.stage pdl1 i.practice_type diag_year age_at_diagnosis  pdl1_given"
 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots

  
 import delimited "all_data_365.csv", clear 
  gen therapy_type = -1
 replace therapy_type = 0 if first_line_chemo == 1
 replace therapy_type = 1 if io_mono == 1
 replace therapy_type = 2 if combo_therapy == 1
 replace therapy_type = 3 if secondary_chemo_drug == 1
 replace therapy_type = 4 if alk_drug == 1
 replace therapy_type = 5 if egfr_drug == 1
 replace therapy_type = 6 if braf_drug == 1
 replace therapy_type = 7 if ros1_drug == 1
 replace therapy_type = 8 if ras_drug == 1
 replace therapy_type = 9 if other_first_line_therapy == 1

   gen insured = 0
 replace insured = 1 if pt_assistance == 1 | other_gov_insurance == 1 | medicare == 1 | medicaid == 1 | commercial_health_plan ==1

 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if therapy_type>=2
 
 ivregress 2sls progression_12 (therapy_type = i.practiceid ) $indiv_covar insured, first robust
 ivreghdfe progression_12  (therapy_type = i.practiceid) $indiv_covar insured, first robust
 
 
 
 ////// Regression discontinuity ///////////
import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/all_data_365.csv", clear
gen therapy_type = -1
replace therapy_type = 0 if first_line_chemo == 1
replace therapy_type = 1 if io_mono == 1
drop if io_mono_used == 1 | io_mono_used>=5
gen first_line = 0
replace first_line = 1 if therapy_type ==1
keep if pdl1_given==1 & therapy_type>=0
  gen insured = 0
 replace insured = 1 if pt_assistance == 1 | other_gov_insurance == 1 | medicare == 1 | medicaid == 1 | commercial_health_plan ==1
rdplot first_line pdl1, c(0.5)
binscatter first_line pdl1, rd(0.5) yti("Probability of IO Monotherapy Treatment") xti("Recentered PD-L1") 
	
graph export "discontinuity_treat.png", replace 

graph twoway (hist pdl1) , xline(0.5, lcolor(red))

kdensity pdl1 , xline(0.5)

//Plotting all, testing only within the optimal bandwidth estimated
rddensity pdl1 , pl c(0.5)

/* Actually do statistical analysis for RD without nivolumab for IO vs chemo */

import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/all_data_365.csv", clear
gen therapy_type = -1
replace therapy_type = 0 if first_line_chemo == 1
replace therapy_type = 1 if io_mono == 1
drop if io_mono_used == 1 | io_mono_used>=5
keep if pdl1_given==1 & therapy_type>=0
gen first_line = 0
replace first_line = 1 if therapy_type ==1
gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5
gen insured = 0
replace insured = 1 if pt_assistance == 1 | other_gov_insurance == 1 | medicare == 1 | medicaid == 1 | commercial_health_plan ==1
drop if alk==1
drop if egfr==1
drop if ros1==1
rdrandinf progression_12 pdl1, cutoff(0.5) fuzzy(first_line itt) covariates(practice_type diag_year age_at_diagnosis insured ecog histology stage race) wl(0.4) wr(0.6) seed(0)
rdrandinf progression_12 pdl1, cutoff(0.5) fuzzy(first_line tsls) covariates(practice_type diag_year age_at_diagnosis insured ecog histology stage race) wl(0.4) wr(0.6) seed(0) firststage
rddensity pdl1, c(0.5) // check sorting/bunching assumption





  