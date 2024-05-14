   **************************** 
 /* Propensity score matching with all DAG factors --> but this time, kaplan-meier curves show days from treatment initiation until progression*/
 
 /*wiresidence mnresidence inresidence varesidence prresidence dcresidence utresidence idresidence moresidence ctresidence nhresidence caresidence arresidence nvresidence deresidence mdresidence tnresidence alresidence njresidence paresidence nyresidence neresidence waresidence wvresidence waresidence azresidence laresidence orresidence okresidence txresidence coresidence iaresidence msresidence riresidence ohresidence scresidence garesidence miresidence ncresidence meresidence flresidence ilresidence nmresidence hiresidence ksresidence kyresidence maresidence */
 ****************************

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
 
  global indiv_covar_no_state "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
 

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
 count if stageii==1 | stageiia ==1 | stageiib ==1  
 count if stageiii ==1 | stageiiia ==1 | stageiiib | stageiiic ==1
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
 graph export "PD_L1_distribution_whole_dataset_include_0_prog.png", replace
 
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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 logit endpoint i.therapy_type ${indiv_covar} i.therapy_type#c.pdl1 alk egfr  ros1  bevacizumabused threeormorechemotherapydrugs  
 predict logit_pred
 
  logit endpoint i.therapy_type ${indiv_covar} i.therapy_type#c.pdl1 alk egfr braf  ros1 kras   bevacizumabused threeormorechemotherapydrugs  if pdl1reported==1
    logit endpoint i.therapy_type ${indiv_covar} i.therapy_type#over_threshold alk egfr braf  ros1 kras   bevacizumabused threeormorechemotherapydrugs  if pdl1reported==1

 
 
  stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Platinum-Based Chemotherapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Other Chemotherapy" 5 "ALK Drug" 6 "EGFR Drug" 7 "BRAF Drug" 8 "ROS1 Drug" 9  "Other First-Line Therapy" 10 "TRK Inhibitor" 11 "MET Drug" 12 "Carboplatin Only" 13 "Cisplatin Only"))
 graph export "prog_survival_with_mutations_pre_match_prog.png", replace
 sts test therapy_type, logrank
 
 drop if therapy_type>=2
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 
 logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} braf  kras  if pdl1reported==1
 logit endpoint i.therapy_type  i.therapy_type#i.over_threshold ${indiv_covar}  braf   kras  if pdl1reported==1
 logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} braf  kras  if pdl1>0
 logit endpoint i.therapy_type  i.therapy_type#i.over_threshold ${indiv_covar} braf   kras  if pdl1>0

 
 
 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation  (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "prog_survival_without_mutations_pre_match_no_combo_prog.png", replace
 sts test therapy_type, logrank

 /*
  drop if hispanicrace==1 | stageia1==1 | connectivetissuedisease==1 | scleroderma == 1| othercnsmetastases == 1 | mnresidence == 1 | dcresidence == 1 | coresidence == 1 | riresidence == 1 | kyresidence == 1 | maresidence == 1 | stage0==1 | occult == 1 | priorkidneytransplant == 1 | priorlivertransplant == 1 | wiresidence == 1 | wvresidence== 1 | orresidence == 1 | deresidence == 1
 teffects psmatch (endpoint)  (therapy_type $indiv_covar pdl1reported)
teffects ipwra (endpoint $indiv_covar pdl1reported)  (therapy_type $indiv_covar pdl1reported), osample(new_var)
teffects ipwra (endpoint $indiv_covar pdl1reported)  (therapy_type $indiv_covar pdl1reported)
teffects ipwra (endpoint $indiv_covar)  (therapy_type $indiv_covar) if pdl1>=0.5
teffects psmatch (endpoint) (therapy_type $indiv_covar) if pdl1>=0.5
teffects ipwra (endpoint $indiv_covar)  (therapy_type $indiv_covar) if pdl1>=0.01 & pdl1<0.5
drop if stageia3==1 | ctresidence == 1 | ksresidence == 1
teffects psmatch (endpoint) (therapy_type $indiv_covar) if pdl1>=0.01 & pdl1<0.5

teffects ipwra (endpoint $indiv_covar)  (therapy_type $indiv_covar) if pdl1==0.0 & pdl1reported==1
*/
logit therapy_type ${indiv_covar}   pdl1reported
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

psmatch2 therapy_type  ${indiv_covar} pdl1reported , outcome(endpoint) caliper(0.2) n(1) noreplacement 
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
  
  
  
     **************************** 
 /* Propensity score matching with all DAG factors --> but this time, kaplan-meier curves show days from treatment initiation until progression*/
 ****************************


/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
 

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
 count if stageii==1 | stageiia ==1 | stageiib ==1  
 count if stageiii ==1 | stageiiia ==1 | stageiiib | stageiiic ==1
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

 
 
  gen time_limit = 365
 gen outcome = "mortality"

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
 graph export "PD_L1_distribution_whole_dataset_include_0_overall.png", replace
 
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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 logit endpoint i.therapy_type ${indiv_covar} i.therapy_type#c.pdl1 alk egfr  ros1  bevacizumabused threeormorechemotherapydrugs  
 predict logit_pred
 //rocreg endpoint logit_pred
 
  logit endpoint i.therapy_type ${indiv_covar} i.therapy_type#c.pdl1 alk egfr braf  ros1 kras   bevacizumabused threeormorechemotherapydrugs  if pdl1reported==1
    logit endpoint i.therapy_type ${indiv_covar} i.therapy_type#over_threshold alk egfr braf  ros1 kras   bevacizumabused threeormorechemotherapydrugs  if pdl1reported==1

 
 
  stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Platinum-Based Chemotherapy" 2 "First-Line IO Monotherapy" 3 "Combination Therapy" 4 "Other Chemotherapy" 5 "ALK Drug" 6 "EGFR Drug" 7 "BRAF Drug" 8 "ROS1 Drug" 9  "Other First-Line Therapy" 10 "TRK Inhibitor" 11 "MET Drug" 12 "Carboplatin Only" 13 "Cisplatin Only"))
 graph export "prog_survival_with_mutations_pre_match_overall.png", replace
 sts test therapy_type, logrank
 
 drop if therapy_type>=2
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 
 logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras  if pdl1reported==1
 logit endpoint i.therapy_type  i.therapy_type#i.over_threshold ${indiv_covar} alk egfr braf  ros1 kras  if pdl1reported==1
 logit endpoint i.therapy_type  i.therapy_type#c.pdl1 ${indiv_covar} alk egfr braf  ros1 kras  if pdl1>0
 logit endpoint i.therapy_type  i.therapy_type#i.over_threshold ${indiv_covar} alk egfr braf  ros1 kras  if pdl1>0

 
 
 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for All Patients in Dataset") subtitle("by Therapy, Pre-Matching") xtitle ("Survival Time from Treatment Initiation  (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "prog_survival_without_mutations_pre_match_no_combo_overall.png", replace
 sts test therapy_type, logrank

 /*
 drop if hispanicrace==1 | stageia1==1 | connectivetissuedisease==1 | scleroderma == 1| othercnsmetastases == 1 | mnresidence == 1 | dcresidence == 1 | coresidence == 1 | riresidence == 1 | kyresidence == 1 | maresidence == 1 | stage0==1 | occult == 1 | priorkidneytransplant == 1 | priorlivertransplant == 1 | wiresidence == 1 | wvresidence== 1 | orresidence == 1 | deresidence == 1
 teffects ipwra (endpoint $indiv_covar pdl1reported)  (therapy_type $indiv_covar pdl1reported)
teffects ipwra (endpoint $indiv_covar)  (therapy_type $indiv_covar) if pdl1>=0.5
teffects psmatch (endpoint) (therapy_type $indiv_covar) if pdl1>=0.5
teffects ipwra (endpoint $indiv_covar)  (therapy_type $indiv_covar) if pdl1>=0.01 & pdl1<0.5
teffects ipwra (endpoint $indiv_covar)  (therapy_type $indiv_covar) if pdl1==0.0 & pdl1reported==1
*/
logit therapy_type ${indiv_covar}  pdl1reported
predict yhat

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))
graph export "propensity_score_pre_matching_overall.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))
graph export "propensity_pre_match_hist_overall.png", replace

pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} pdl1reported , outcome(endpoint) caliper(0.2) n(1) noreplacement 
count 


sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))
graph export "propensity_post_match_overall.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 
graph export "propensity_post_match_hist_overall.png", replace



 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for Entire Dataset") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_overall.png", replace
  sts test therapy_type, logrank
  
  
  /////////////
  

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused  braf kras"
 

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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 
 
 drop if therapy_type>=2
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 keep if pdl1>=0.5


logit therapy_type ${indiv_covar}  pdl1reported
predict yhat

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))


pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} pdl1reported , outcome(endpoint) caliper(0.2) n(1) noreplacement 
count 


sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 

 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset (PD-L1 >= 50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_prog_only_match_on_over50.png", replace
  sts test therapy_type, logrank
  

///////////////////

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused  braf kras"
 

 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots
import delimited "all_data_365_10000.csv", clear 
 gen time_limit = 365
 gen outcome = "mortality"
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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 
 
 drop if therapy_type>=2
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 keep if pdl1>=0.5


logit therapy_type ${indiv_covar}  pdl1reported
predict yhat

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))


pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} pdl1reported , outcome(endpoint) caliper(0.2) n(1) noreplacement 
count 


sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 

 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for Entire Dataset (PD-L1 >= 50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_overall_only_match_on_over50.png", replace
  sts test therapy_type, logrank
  


  ///////////////////

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
 

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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 
 
 drop if therapy_type>=2
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 keep if pdl1==0.0 & pdl1reported==1


logit therapy_type ${indiv_covar}  pdl1reported
predict yhat

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))


pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} pdl1reported , outcome(endpoint) caliper(0.2) n(1) noreplacement 
count 


sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 

 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset (PD-L1 Negative)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_prog_only_match_on_confirmed_negative.png", replace
  sts test therapy_type, logrank
  

  

  ///////////////////

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
 

 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots
import delimited "all_data_365_10000.csv", clear 
 gen time_limit = 365
 gen outcome = "mortality"
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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 
 
 drop if therapy_type>=2
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 keep if pdl1==0.0 & pdl1reported==1


logit therapy_type ${indiv_covar}  pdl1reported
predict yhat

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))


pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} pdl1reported , outcome(endpoint) caliper(0.2) n(1) noreplacement 
count 


sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 

 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for Entire Dataset (PD-L1 Negative)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_overall_only_match_on_confirmed_negative.png", replace
  sts test therapy_type, logrank
  
  
  
  
  
  
  
  
    /////////////////// 1<PDL1<50

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
 

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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 
 
 drop if therapy_type>=2
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 keep if pdl1>=0.01 & pdl1<0.5


logit therapy_type ${indiv_covar}  pdl1reported
predict yhat

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))


pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} pdl1reported , outcome(endpoint) caliper(0.2) n(1) noreplacement 
count 


sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 

 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival for Entire Dataset (1% <= PD-L1 < 50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_prog_only_match_on_1_50.png", replace
  sts test therapy_type, logrank
  

  

  ///////////////////

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
 

 global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
 cd "${path}"
 set scheme cleanplots
import delimited "all_data_365_10000.csv", clear 
 gen time_limit = 365
 gen outcome = "mortality"
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

 
 replace censor_time  = time_limit if censor_time == 0 | censor_time > time_limit
 
 
 drop if therapy_type>=2
 drop if alk==1
 drop if egfr==1
 drop if ros1==1
 drop if bevacizumabused == 1 | threeormorechemotherapydrugs == 1
 drop if clinicalstudydrugused == 1
 keep if pdl1>=0.01 & pdl1<0.5


logit therapy_type ${indiv_covar}  pdl1reported
predict yhat

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "First-Line IO Monotherapy"))


pstest ${indiv_covar}, raw treated(therapy_type)

psmatch2 therapy_type  ${indiv_covar} pdl1reported , outcome(endpoint) caliper(0.2) n(1) noreplacement 
count 


sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${indiv_covar} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemotherapy") label(2 "IO Monotherapy"))


graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "First-Line Chemotherapy") ///
		label(2 "IO Monotherapy")) 

 stset censor_time, failure(endpoint)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival for Entire Dataset (1% <= PD-L1 < 50%)") subtitle("by Therapy (Post-Matching)") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "First-Line Chemotherapy" 2 "First-Line IO Monotherapy")) 
 graph export "chemo_io_post_match_overall_only_match_on_1_50.png", replace
  sts test therapy_type, logrank
  
