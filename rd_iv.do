
 ////// Regression discontinuity ///////////

  global indiv_covar_no_pdl1 "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
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

drop if glucocorticoidusepriortotreatmen==1
drop if threeormorechemotherapydrugs == 1
drop if bevacizumabused==1
keep if pdl1 > 0

 gen insured = 0
 replace insured = 1 if patientassistanceprogram == 1 | othergovernmentalinsurance == 1 | medicare == 1 | medicaid == 1 | commercialhealthplan ==1
replace daysfromadvanceddiagnosistotreat = log(daysfromadvanceddiagnosistotreat)
replace ageatdiagnosis = log(ageatdiagnosis)
/*
 binscatter squamouscellcarcinoma pdl1, yti("Proportion with Squamous Cell Carcinoma") xti("PD-L1") 
graph export "squamous_pdl1.png", replace */

rdplot first_line pdl1, c(0.5) 
graph export "polynomial_fit_RD.png", replace
binscatter first_line pdl1, rd(0.5) yti("Proportion Receiving IO Monotherapy Treatment") xti("PD-L1") 
graph export "discontinuity_treat.png", replace 
/*
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
 */

kdensity pdl1 , xline(0.5)

//Plotting all, testing only within the optimal bandwidth estimated
rddensity pdl1 , pl c(0.5)


// Note: The "T" is your local average treatment effect. The P>|T| in the large sample is your p-value. Use the confidence interval.
// RD density plot shows that while there is a significant spike from 0.4 to 0.5 there is no reason to think pathologists
// are artificially inflating PD-L1 values so as to increase likelihood of IO monotherapy. 
/* Actually do statistical analysis for RD without nivolumab for IO vs chemo */
rddensity pdl1, c(0.5) vce(jackknife) plot

rdplot progression_outcome  pdl1, c(0.5) 
rdwinselect pdl1 $indiv_covar_no_pdl1, c(0.5) seed(0) reps(1000) level(0.05) wmass

rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.4) wr(0.51)  firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 


drop if therapyyear < 2014 

rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 


  /////////////////////////////
/*  Instrumental variables */

global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
cd "${path}"
set scheme cleanplots
 global indiv_covar "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma pdl1 hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused  braf kras"
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

