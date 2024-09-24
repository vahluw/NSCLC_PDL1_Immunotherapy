/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */



global therapy_type_decider "sdh1 sdh2 sdh3 sdh4 sdh5 braf kras ecog0 ecog1 ecog2 ecog3 ecog4 nonsquamouscellcarcinoma squamouscellcarcinoma pdl1  hispanicethnicity  diagnosisyear ageatdiagnosis white asian black otherrace hispanicrace male female  daysfromadvanceddiagnosistotreat patientassistanceprogram othergovernmentalinsurance medicare selfpay medicaid commercialhealthplan noinsurance  neversmoker  academicmedicalcenter kidney_bool liver_bool connective_tissue_bool diabetes bonemetastases brainmetastases digestivesystemmetastases adrenalmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen over_threshold"


global therapy_type_decider_categorical "i.sdh braf kras i.ecog i.histology i.race  hispanicethnicity  diagnosisyear age_dichotomized i.gender    medicare  medicaid commercialhealthplan otherinsurance noinsurance  i.smokinghistory  academicmedicalcenter   bonemetastases brainmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "

global therapy_type_no_mut "i.sdh i.ecog i.histology i.race  hispanicethnicity  diagnosisyear age_dichotomized i.gender    medicare  medicaid commercialhealthplan otherinsurance noinsurance  i.smokinghistory  academicmedicalcenter   bonemetastases brainmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "

global therapy_type_no_pdl1 "i.sdh braf kras i.ecog i.histology i.race   hispanicethnicity  diagnosisyear ageatdiagnosis i.gender  daysfromadvanceddiagnosistotreat   medicare selfpay medicaid commercialhealthplan otherinsurance uninsured  i.smokinghistory  academicmedicalcenter kidney_bool liver_bool connective_tissue_bool diabetes bonemetastases brainmetastases digestivesystemmetastases adrenalmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "

global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/research_letter/"
cd "${path}"
set scheme cleanplots

import delimited "all_data.csv", clear 

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
gen kidney_bool = (chronickidneydisease == 1 |  priorkidneytransplant == 1)
gen liver_bool = (cirrhosis == 1 | hepatitis == 1 | priorlivertransplant  == 1)
replace otherfirstlinetherapy = 1 if antirasdrug==1

gen ecog = 2
replace ecog = 0 if ecog0==1
replace ecog = 1 if ecog1==1


tab connective_tissue_bool
tab kidney_bool
tab liver_bool

gen time_limit = 365
gen outcome = "progression"

//drop if ecog0==0 & ecog1 == 0 & ecog2 == 0 & ecog3==0 & ecog4 == 0

replace daysfromadvanceddiagnosistotreat = log(daysfromadvanceddiagnosistotreat)
gen insured = 0
replace insured = 1 if patientassistanceprogram == 1 | othergovernmentalinsurance == 1 | medicare == 1 | medicaid == 1 | commercialhealthplan ==1


//replace ageatdiagnosis = log(ageatdiagnosis)

gen therapy_type = 1
replace therapy_type = 0 if firstlinechemotherapy == 1
replace therapy_type = 2 if firstlinecombinationtherapy == 1
replace therapy_type = 3 if nonfirstlinechemotherapy == 1
replace therapy_type = 4 if antialkdrug == 1
replace therapy_type = 5 if antiegfrdrug == 1
replace therapy_type = 6 if antibrafdrug == 1
replace therapy_type = 7 if antiros1drug == 1
replace therapy_type = 8  if otherfirstlinetherapy == 1 
replace therapy_type = 9 if trkinhibitor == 1
replace therapy_type = 10 if metinhibitor== 1
replace therapy_type = 11 if carboplatinmonotherapy == 1
replace therapy_type = 12 if cisplatinmonotherapy == 1

gen io_mono = (therapy_type==1)

drop if therapy_type > 2

drop if pdl1reported==0

drop if alk==1 | egfr == 1 | ros1 == 1 | therapy_type > 2

gen sdh1 = 0
gen sdh2 = 0
gen sdh3 = 0
gen sdh4 = 0
gen sdh5 = 0
replace sdh1 = 1 if sdh==1
replace sdh2 = 1 if sdh==2
replace sdh3 = 1 if sdh==3
replace sdh4 = 1 if sdh==4
replace sdh5 = 1 if sdh==5
drop if therapy_type==0
replace therapy_type = 0 if therapy_type==2

gen histology = -1
replace histology = 0 if squamouscellcarcinoma == 1
replace histology = 1 if nonsquamouscellcarcinoma == 1

binscatter therapy_type  pdl1 if therapy_type  <=1, ytitle ("Proportion Receiving IO Monotherapy") xtitle("PD-L1 Expression")
tab therapy_type diagnosisyear
reg therapy_type  pdl1 if therapy_type  <=1

drop if ecog0 == 0 & ecog1==0 & ecog2==0 & ecog3 == 0 & ecog4==0
drop if white == 0 & black == 0 & hispanicrace == 0 & asian == 0 & otherrace==0
drop if sdh==0
drop if nonsquamouscellcarcinoma== 0 & squamouscellcarcinoma==0
drop if male == 0 & female ==0
drop if neversmoker==0 & previoussmoker==0

gen race = -1
replace race = 0 if white == 1
replace race = 1 if black == 1
replace race = 2 if asian == 1
replace race = 3 if hispanicrace == 1
replace race = 4 if otherrace == 1

drop if creatinine == 0.0
drop if ast == 0 & alt == 0 & bilirubin == 0.0

replace estimated_gfr = log(estimated_gfr)
replace ast = log(ast)
replace alt = log(alt)
replace bilirubin = log(bilirubin)

gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5

gen gender = -1
replace gender = 0 if male == 1
replace gender = 1 if female == 1

keep if pdl1reported==1

gen smokinghistory = -1
replace smokinghistory = 0 if neversmoker == 1
replace smokinghistory = 1 if previoussmoker==1

gen uninsured = 0
replace uninsured = 1 if (noinsurance == 1 | selfpay==1) & medicare == 0 & medicaid == 0 & commercialhealthplan==0 & patientassistanceprogram== 0 & othergovernmentalinsurance== 0

gen otherinsurance = 0
replace otherinsurance = 1 if  (patientassistanceprogram== 1 | othergovernmentalinsurance== 1)

drop if pdl1 < 0.01
drop if diagnosisyear < 2018
gen age_dichotomized = 0
replace age_dichotomized = 1 if ageatdiagnosis>=75
binscatter therapy_type  diagnosisyear if therapy_type  <=1, ytitle ("Proportion Receiving IO Monotherapy") xtitle("Year of Advanced Diagnosis")
binscatter therapy_type  diagnosisyear if therapy_type  <=1 & pdl1>=0.5, ytitle ("Proportion Receiving IO Monotherapy") xtitle("Year of Advanced Diagnosis")
binscatter therapy_type  diagnosisyear if therapy_type  <=1, by (over_threshold) ytitle ("Proportion Receiving IO Monotherapy") xtitle("Year of Advanced Diagnosis") legend(order(1 "1% <= PD-L1 < 50%" 2 "PD-L1 >= 50%"))


logistic therapy_type  $therapy_type_decider_categorical pdl1
logistic therapy_type  $therapy_type_decider_categorical if pdl1>=0.5
logistic therapy_type  $therapy_type_decider_categorical if pdl1<0.5 

logistic therapy_type  $therapy_type_no_mut if pdl1>=0.5
logistic therapy_type  $therapy_type_no_mut if pdl1<0.5 

count if pdl1>=0.5
sum ageatdiagnosis if pdl1>=0.5
tab sdh if pdl1>=0.5
tab gender if pdl1>=0.5
tab race if pdl1>=0.5
tab hispanicethnicity if pdl1>=0.5
tab ecog if pdl1>=0.5
tab academicmedicalcenter if pdl1>=0.5
tab smokinghistory if pdl1>=0.5
tab histology if pdl1>=0.5
tab medicare if pdl1>=0.5
tab medicaid if pdl1>=0.5
tab commercialhealthplan if pdl1>=0.5
tab otherinsurance if pdl1>=0.5
tab noinsurance if pdl1>=0.5
tab therapy_type if pdl1>=0.5
tab braf if pdl1>=0.5
tab kras if pdl1>=0.5
tab antiinfectiveusepriortotreatment if pdl1>=0.5
tab glucocorticoidusepriortotreatmen if pdl1>=0.5



count if pdl1 < 0.5
sum ageatdiagnosis if pdl1 < 0.5
tab sdh if pdl1 < 0.5
tab gender if pdl1 < 0.5
tab race if pdl1 < 0.5
tab hispanicethnicity if pdl1 < 0.5
tab ecog if pdl1 < 0.5
tab academicmedicalcenter if pdl1 < 0.5
tab smokinghistory if pdl1 < 0.5
tab histology if pdl1 < 0.5
tab medicare if pdl1 < 0.5
tab medicaid if pdl1 < 0.5
tab commercialhealthplan if pdl1 < 0.5
tab otherinsurance if pdl1 < 0.5
tab noinsurance if pdl1 < 0.5
tab therapy_type if pdl1 < 0.5
tab braf if pdl1 < 0.5
tab kras if pdl1 < 0.5
tab antiinfectiveusepriortotreatment if pdl1 < 0.5
tab glucocorticoidusepriortotreatmen if pdl1 < 0.5
