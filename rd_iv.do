
 ////// Regression discontinuity ///////////

  global indiv_covar_no_pdl1 "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
  
  global psm "braf kras ecog0 ecog1 ecog2 ecog3 ecog4 nonsquamouscellcarcinoma  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance academicmedicalcenter kidney_bool liver_bool connective_tissue_bool interstitiallungdisease  bonemetastases brainmetastases digestivesystemmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen"
  
  
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
 gen therapy_type = 1
 replace therapy_type = 2 if firstlinechemotherapy == 1
 replace therapy_type = 0 if firstlinecombinationtherapy == 1
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
rdwinselect pdl1 $psm if therapyyear>=2019, c(0.5) seed(0) reps(1000) level(0.05) wmass

rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.3) wr(0.6)  firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 


drop if therapyyear < 2014 

rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 


//////





 ////// Regression discontinuity ///////////

  global indiv_covar_no_pdl1 "ecog0 ecog1 ecog2 ecog3 ecog4 squamouscellcarcinoma nonsquamouscellcarcinoma hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance stage0 stageia stageia1 stageia2 stageia3 stageib stageii stageiia stageiib stageiii stageiiia stageiiib stageiiic stageiv stageiva stageivb occult neversmoker  academicmedicalcenter chronickidneydisease  priorkidneytransplant cirrhosis hepatitis priorlivertransplant connectivetissue scleroderma lupus rheumatoidarthritis interstitiallungdisease diabetes bonemetastases brainmetastases othercnsmetastases digestivesystemmetastases adrenalmetastases unspecifiedmetastases  clinicalstudydrug creatinine bilirubin ast alt albumin antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen clinicalstudydrugused braf kras"
  
  global psm "braf kras ecog0 ecog1 ecog2 ecog3 ecog4 nonsquamouscellcarcinoma  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance academicmedicalcenter kidney_bool liver_bool connective_tissue_bool interstitiallungdisease  bonemetastases brainmetastases digestivesystemmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen"
  
  
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
rdwinselect pdl1 $psm if nonsquamouscellcarcinoma==1 & therapyyear>=2016, c(0.5) seed(0) reps(1000) level(0.05) wmass

rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.01) wr(1)  firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 


drop if therapyyear < 2014 

rdrandinf progression_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 
rdrandinf mortality_outcome pdl1, cutoff(0.5) fuzzy(first_line tsls) kernel(uniform) seed(0)  ci(.05) wl (0.1) wr(0.8)  firststage 

  /////////////////////////////
/*  Instrumental variables */

global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"
cd "${path}"
set scheme cleanplots
 global indiv_covar "braf kras ecog0 ecog1 ecog2 ecog3 ecog4 nonsquamouscellcarcinoma  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance academicmedicalcenter kidney_bool liver_bool connective_tissue_bool interstitiallungdisease  bonemetastases brainmetastases digestivesystemmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen"
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
 gen therapy_type = 1
 replace therapy_type = 2 if firstlinechemotherapy == 1
 replace therapy_type = 0 if firstlinecombinationtherapy == 1
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
keep if therapyyear>=2019 & pdl1>0
gen over_threshold =(pdl1>=0.5)
//keep if pdl1_given==1
ivreg2 progression_outcome (therapy_type=i.practiceid) // Unadjusted regression
ivreg2 progression_outcome (therapy_type=i.practiceid ) $indiv_covar   // Adjusted regression

keep if pdl1reported==1
ivreg2 progression_outcome (therapy_type=over_threshold) // Unadjusted regression
ivreg2 progression_outcome (therapy_type=over_threshold i.practiceid ) $indiv_covar // Adjusted regression

