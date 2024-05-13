   **************************** 
 /* Propensity score matching with all DAG factors --> but this time, kaplan-meier curves show days from treatment initiation until progression*/
 ****************************

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male female daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker communitymedicalcenter academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused"

 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots
 



 import delimited "all_data_365_10000.csv", clear 

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
 graph export "PD_L1_distribution_whole_dataset_include_0.png", replace
 
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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 logit endpoint i.therapy_type ${indiv_covar} alk egfr braf  ros1 kras   bevacizumabused threeormorechemotherapydrugs  
 //rocreg endpoint logit_pred
 
 
  stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Platinum-Based Chemotherapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Other Chemotherapy" 5 "ALK Drug" 6 "EGFR Drug" 7 "BRAF Drug" 8 "ROS1 Drug" 9  "Other First-Line Therapy" 10 "TRK Inhibitor" 11 "MET Drug" 12 "Carboplatin Only" 13 "Cisplatin Only"))
 graph export "prog_survival_with_mutations_pre_match_prog.png", replace
 sts test therapy_type, logrank
 
  logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras  if pdl1reported==1
 drop if therapy_type>=2

  logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras  if pdl1reported==1


 drop if alk==1
 drop if egfr==1
 drop if ros1==1

 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 
   logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras   if pdl1reported==1
   
      logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras  if pdl1 > 0

 
 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation  (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "prog_survival_without_mutations_pre_match_no_combo_prog.png", replace
 sts test therapy_type, logrank


 logit therapy_type ${indiv_covar} braf kras   //i.therapy_type##over_threshold
 
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

psmatch2 therapy_type  ${indiv_covar} pdl1reported braf kras , outcome(endpoint) caliper(0.2) //n(1) noreplacement 
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
  
     stset censor_time if pdl1==0.0 & pdl1reported==1, failure(endpoint)
 stci if pdl1==0.0 & pdl1reported==1, by(therapy_type) rmean
 stcox therapy_type if pdl1==0.0 & pdl1reported==1
 sts graph if pdl1==0.0 & pdl1reported==1, by(therapy_type) title("Progression-Free Survival for Entire Dataset (Confirmed PD-L1 Negative)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
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
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male female daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker communitymedicalcenter academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused"

 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots


 import delimited "all_data_365_10000.csv", clear 

  gen time_limit = 365
 gen outcome = "mortality"

 
 
 gen censor_time = progression_days 
 replace censor_time = mortality_days if outcome == "mortality"
 gen endpoint = progression_outcome
 replace endpoint = mortality_outcome if outcome == "mortality"
 


 tab progression_outcome
 gen insured = 0
 replace insured = 1 if patientassistanceprogram == 1 | othergovernmentalinsurance == 1 | medicare == 1 | medicaid == 1 | commercialhealthplan ==1


 histogram pdl1, percent bin(9) xtitle("PD-L1 Intensity") color(ebblue)
 graph export "PD_L1_distribution_whole_dataset_include_0.png", replace
 
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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 logit endpoint i.therapy_type ${indiv_covar} alk egfr braf  ros1 kras   bevacizumabused threeormorechemotherapydrugs  
 //rocreg endpoint logit_pred
 
 
  stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Platinum-Based Chemotherapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Other Chemotherapy" 5 "ALK Drug" 6 "EGFR Drug" 7 "BRAF Drug" 8 "ROS1 Drug" 9  "Other First-Line Therapy" 10 "TRK Inhibitor" 11 "MET Drug" 12 "Carboplatin Only" 13 "Cisplatin Only"))
 graph export "prog_survival_with_mutations_pre_match_prog.png", replace
 sts test therapy_type, logrank
 
  logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras  if pdl1reported==1
 drop if therapy_type>=2

  logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras  if pdl1reported==1


 drop if alk==1
 drop if egfr==1
 drop if ros1==1

 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 
   logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras   if pdl1reported==1
   
      logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras  if pdl1 > 0

 
 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation  (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "prog_survival_without_mutations_pre_match_no_combo_prog.png", replace
 sts test therapy_type, logrank


 logit therapy_type ${indiv_covar} braf kras   //i.therapy_type##over_threshold
 
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

psmatch2 therapy_type  ${indiv_covar} pdl1reported braf kras , outcome(endpoint) caliper(0.2) //n(1) noreplacement 
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
  
     stset censor_time if pdl1==0.0 & pdl1reported==1, failure(endpoint)
 stci if pdl1==0.0 & pdl1reported==1, by(therapy_type) rmean
 stcox therapy_type if pdl1==0.0 & pdl1reported==1
 sts graph if pdl1==0.0 & pdl1reported==1, by(therapy_type) title("Progression-Free Survival for Entire Dataset (Confirmed PD-L1 Negative)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_pdl1_zero_prog.png", replace
  sts test therapy_type, logrank
  
 keep if pdl1>=0.01 & pdl1 < 0.5
 stset censor_time , failure(endpoint)
 stci , by(therapy_type) rmean
 stcox therapy_type 
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset (1%<=PD-L1<50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_pdl1_under50_prog.png", replace
  sts test therapy_type, logrank
  
 ////// Regression discontinuity ///////////

import delimited "all_data_365_10000.csv", clear
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
 keep if therapy_type == 0 | therapy_type==1
 drop if firstlinenivolumabmonotherapy==1
 drop if firstlineipilimumabnivolumab ==1

 drop if alk==1
drop if egfr==1
drop if ros1==1
gen first_line = 0
replace first_line = 1 if therapy_type ==1
gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5

drop if threeormorechemotherapydrugs == 1
drop if bevacizumabused==1
keep if pdl1> 0
keep if pdl1reported==1

 gen insured = 0
 replace insured = 1 if patientassistanceprogram == 1 | othergovernmentalinsurance == 1 | medicare == 1 | medicaid == 1 | commercialhealthplan ==1


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
binscatter insured pdl1, rd(0.5) yti("Proportion Insured") xti("PD-L1") 
graph export "discontinuity_insured.png", replace 
binscatter never_smoker pdl1, rd(0.5) yti("Proportion Never Smoker") xti("PD-L1") 
graph export "discontinuity_never_smoker.png", replace
binscatter prev_smoker pdl1, rd(0.5) yti("Proportion Previous Smoker") xti("PD-L1") 
graph export "discontinuity_prev_smoker.png", replace
 
kdensity pdl1 , xline(0.5)

//Plotting all, testing only within the optimal bandwidth estimated
rddensity pdl1 , pl c(0.5)


// Note: The "T" is your local average treatment effect. The P>|T| in the large sample is your p-value. Use the confidence interval.
// RD density plot shows that while there is a significant spike from 0.4 to 0.5 there is no reason to think pathologists
// are artificially inflating PD-L1 values so as to increase likelihood of IO monotherapy. 
/* Actually do statistical analysis for RD without nivolumab for IO vs chemo */
rddensity pdl1, c(0.5) vce(jackknife) plot

rdplot progression_outcome  pdl1, c(0.5) 
rdwinselect pdl1 male daysfromadvanceddiagnosistotreat therapyyear ageatdiagnosis insured neversmoker  black white asian otherrace  academicmedicalcenter, c(0.5) seed(0) reps(1000) level(0.05) wmass
rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 


drop if therapyyear < 2014 

rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 


  /////////////////////////////
/*  Instrumental variables */

global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
cd "${path}"
set scheme cleanplots
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male female daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker communitymedicalcenter academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused"
import delimited "all_data_365_10000.csv", clear 

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
  keep if therapy_type == 0 | therapy_type==1
drop if alk==1
drop if egfr==1
drop if ros1==1
//drop if bev_used == 1 | three_plus_chemo_drugs == 1
keep if therapy_type >= 0
set emptycells drop

logit progression_outcome i.therapy_type $indiv_covar 
gen over_threshold =(pdl1>=0.5)
//keep if pdl1_given==1
ivreg2 progression_outcome (therapy_type=i.practiceid) // Unadjusted regression
ivreg2 progression_outcome (therapy_type=i.practiceid ) $indiv_covar   // Adjusted regression

keep if pdl1reported==1
ivreg2 progression_outcome (therapy_type=over_threshold) // Unadjusted regression
ivreg2 progression_outcome (therapy_type=over_threshold i.practiceid ) $indiv_covar // Adjusted regression

