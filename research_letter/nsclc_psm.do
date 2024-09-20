/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */



//global therapy_type_no_pdl1 "i.sdh braf kras i.ecog i.histology i.race   hispanicethnicity  diagnosisyear ageatdiagnosis i.gender  daysfromadvanceddiagnosistotreat   medicare selfpay medicaid commercialhealthplan otherinsurance uninsured  i.smokinghistory  academicmedicalcenter diabetes bonemetastases brainmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "

global therapy_type_no_pdl1 "i.sdh braf kras i.ecog i.histology i.race diagnosisyear ageatdiagnosis i.gender  daysfromadvanceddiagnosistotreat  i.smokinghistory  academicmedicalcenter antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "


global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/research_letter/"
cd "${path}"
set scheme cleanplots

import delimited "all_data.csv", clear 

gen ecog = 5
replace ecog = 0 if ecog0==1
replace ecog = 1 if ecog1==1
replace ecog =2 if ecog2 == 1 
replace ecog = 3 if ecog3==1
replace ecog = 4 if ecog4==1

gen time_limit = 365
gen outcome = "progression"


gen insured = 0
replace insured = 1 if patientassistanceprogram == 1 | othergovernmentalinsurance == 1 | medicare == 1 | medicaid == 1 | commercialhealthplan ==1

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

gen histology = 2
replace histology = 0 if squamouscellcarcinoma == 1
replace histology = 1 if nonsquamouscellcarcinoma == 1

gen race = 5
replace race = 0 if white == 1
replace race = 1 if black == 1
replace race = 2 if asian == 1
replace race = 3 if hispanicrace == 1
replace race = 4 if otherrace == 1

gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5

gen gender = 2
replace gender = 0 if male == 1
replace gender = 1 if female == 1

keep if pdl1reported==1

gen smokinghistory = 2
replace smokinghistory = 0 if neversmoker == 1
replace smokinghistory = 1 if previoussmoker==1

gen uninsured = 0
replace uninsured = 1 if (noinsurance == 1 | selfpay==1) & medicare == 0 & medicaid == 0 & commercialhealthplan==0 & patientassistanceprogram== 0 & othergovernmentalinsurance== 0

gen otherinsurance = 0
replace otherinsurance = 1 if  (patientassistanceprogram== 1 | othergovernmentalinsurance== 1)

drop if timetocensor < 365 & (progressionoutcome == 0 | mortalityoutcome==0)

 replace progressionoutcome = 1 if progressiondays <= 365  & progressiondays>0
 replace mortalityoutcome = 1 if mortalitydays > 0 & mortalitydays <=365
 replace progressiondays = timetocensor if timetocensor < 365 & progressionoutcome==0
 replace mortalitydays = timetocensor if timetocensor < 365 & mortalityoutcome==0
 replace progressiondays = 365 if progressiondays > 365
 replace mortalitydays = 365 if mortalitydays > 365
   replace progressiondays = 365 if progressiondays == 0 | progressiondays > 365
  replace mortalitydays = 365 if mortalitydays == 0 | mortalitydays > 365

  drop if race == 3 | ecog==4
teffects psmatch (progressionoutcome) (therapy_type $therapy_type_no_pdl1 pdl1)  if pdl1>0,  nneighbor(1)
teffects psmatch (mortalityoutcome) (therapy_type $therapy_type_no_pdl1 pdl1) if pdl1>0,  nneighbor(1)
teffects ipw (progressionoutcome) (therapy_type $therapy_type_no_pdl1 pdl1)  if pdl1>0 
teffects ipw (mortalityoutcome) (therapy_type $therapy_type_no_pdl1 pdl1) if pdl1>0 
 
teffects psmatch (progressionoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 <0.5 & pdl1>0, nneighbor(1)
teffects psmatch (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 <0.5 & pdl1>0, nneighbor(1)
teffects ipw (progressionoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 <0.5 & pdl1>0
teffects ipw (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 <0.5 & pdl1>0
teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1)
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1)

teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1 <0.5 & pdl1>0
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1 <0.5 & pdl1>0

 drop if pdl1>=0.5 | pdl1==0
 logit therapy_type $therapy_type_no_pdl1 if therapy_type <=1
 predict yhat
 

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))
graph export "propensity_score_pre_matching_progression_reck.png", replace
graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "IO Monotherapy") ///
		label(2 "Chemoimmunotherapy"))
graph export "propensity_pre_match_hist_prog_reck.png", replace


 stset progressiondays, failure(progressionoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival") subtitle("1% < PD-L1 < 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_under_50_pfs.png", replace
  sts test therapy_type, logrank
  


 stset mortalitydays, failure(mortalityoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival") subtitle("1% < PD-L1 < 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_under_50_os.png", replace
  sts test therapy_type, logrank
  


pstest ${therapy_type_no_pdl1}, raw treated(therapy_type)

psmatch2 therapy_type  ${therapy_type_no_pdl1} , outcome(progressionoutcome) caliper(0.2) n(1) noreplacement 
count 

sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${therapy_type_no_pdl1} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemoimmunotherapy") label(2 "IO Monotherapy"))
graph export "propensity_post_match_reck.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "Chemoimmunotherapy") ///
		label(2 "IO Monotherapy")) 
graph export "propensity_post_match_hist_prog_over_50.png", replace



 stset progressiondays, failure(progressionoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival") subtitle("1% < PD-L1 < 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_under_50_pfs.png", replace
  sts test therapy_type, logrank
  


 stset mortalitydays, failure(mortalityoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival") subtitle("1% < PD-L1 < 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_under_50_os.png", replace
  sts test therapy_type, logrank
  
  
pstest ${therapy_type_no_pdl1}, treated(therapy_type)

////////

/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */



//global therapy_type_no_pdl1 "i.sdh braf kras i.ecog i.histology i.race   hispanicethnicity  diagnosisyear ageatdiagnosis i.gender  daysfromadvanceddiagnosistotreat   medicare selfpay medicaid commercialhealthplan otherinsurance uninsured  i.smokinghistory  academicmedicalcenter diabetes bonemetastases brainmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "

global therapy_type_no_pdl1 "i.sdh braf kras i.ecog i.histology i.race     diagnosisyear ageatdiagnosis i.gender  daysfromadvanceddiagnosistotreat    i.smokinghistory  academicmedicalcenter  antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "


global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/research_letter/"
cd "${path}"
set scheme cleanplots

import delimited "all_data.csv", clear 

gen ecog = 5
replace ecog = 0 if ecog0==1
replace ecog = 1 if ecog1==1
replace ecog =2 if ecog2 == 1 
replace ecog = 3 if ecog3==1
replace ecog = 4 if ecog4==1


gen time_limit = 365
gen outcome = "progression"


gen insured = 0
replace insured = 1 if patientassistanceprogram == 1 | othergovernmentalinsurance == 1 | medicare == 1 | medicaid == 1 | commercialhealthplan ==1


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

gen histology = 2
replace histology = 0 if squamouscellcarcinoma == 1
replace histology = 1 if nonsquamouscellcarcinoma == 1


drop if timetocensor < 365 & (progressionoutcome == 0 | mortalityoutcome==0)

gen race = 5
replace race = 0 if white == 1
replace race = 1 if black == 1
replace race = 2 if asian == 1
replace race = 3 if hispanicrace == 1
replace race = 4 if otherrace == 1

gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5

gen gender = 2
replace gender = 0 if male == 1
replace gender = 1 if female == 1

keep if pdl1reported==1

gen smokinghistory = 2
replace smokinghistory = 0 if neversmoker == 1
replace smokinghistory = 1 if previoussmoker==1

gen uninsured = 0
replace uninsured = 1 if (noinsurance == 1 | selfpay==1) & medicare == 0 & medicaid == 0 & commercialhealthplan==0 & patientassistanceprogram== 0 & othergovernmentalinsurance== 0

gen otherinsurance = 0
replace otherinsurance = 1 if  (patientassistanceprogram== 1 | othergovernmentalinsurance== 1)


 replace progressionoutcome = 1 if progressiondays <= 365  & progressiondays>0
 replace mortalityoutcome = 1 if mortalitydays > 0 & mortalitydays <=365
 replace progressiondays = timetocensor if timetocensor < 365 & progressionoutcome==0
 replace mortalitydays = timetocensor if timetocensor < 365 & mortalityoutcome==0
 replace progressiondays = 365 if progressiondays > 365
 replace mortalitydays = 365 if mortalitydays > 365
  replace progressiondays = 365 if progressiondays == 0 | progressiondays > 365
  replace mortalitydays = 365 if mortalitydays == 0 | mortalitydays > 365
  
   drop if race == 3 | ecog==4
 teffects psmatch (progressionoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5,  nneighbor(1)
teffects psmatch (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5,  nneighbor(1)
teffects ipw (progressionoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5
teffects ipw (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5
teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5


 drop if pdl1<0.5
 logit therapy_type $therapy_type_no_pdl1 if therapy_type <=1
 predict yhat


graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))
graph export "propensity_score_pre_matching_progression_reck.png", replace
graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "IO Monotherapy") ///
		label(2 "Chemoimmunotherapy"))
graph export "propensity_pre_match_hist_prog_reck.png", replace

 stset progressiondays, failure(progressionoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival") subtitle("PD-L1 >= 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_over_50_pfs.png", replace
  sts test therapy_type, logrank
  


 stset mortalitydays, failure(mortalityoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival") subtitle("PD-L1 >= 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_over_50_os.png", replace
  sts test therapy_type, logrank
  


pstest ${therapy_type_no_pdl1}, raw treated(therapy_type)

psmatch2 therapy_type  ${therapy_type_no_pdl1} , outcome(progressionoutcome) caliper(0.2) n(1) noreplacement 
count 

sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${therapy_type_no_pdl1} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemoimmunotherapy") label(2 "IO Monotherapy"))
graph export "propensity_post_match_reck.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "Chemoimmunotherapy") ///
		label(2 "IO Monotherapy")) 
graph export "propensity_post_match_hist_prog_over_50.png", replace



 stset progressiondays, failure(progressionoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival") subtitle("PD-L1 >= 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_over_50_pfs.png", replace
  sts test therapy_type, logrank
  


 stset mortalitydays, failure(mortalityoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival") subtitle("PD-L1 >= 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_over_50_os.png", replace
  sts test therapy_type, logrank
  
  
pstest ${therapy_type_no_pdl1}, treated(therapy_type)

/////

global therapy_type_no_pdl1 "i.sdh braf kras i.ecog i.histology i.race     diagnosisyear ageatdiagnosis i.gender  daysfromadvanceddiagnosistotreat    i.smokinghistory  academicmedicalcenter  antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "


global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/research_letter/"
cd "${path}"
set scheme cleanplots

import delimited "all_data.csv", clear 

gen ecog = 5
replace ecog = 0 if ecog0==1
replace ecog = 1 if ecog1==1
replace ecog =2 if ecog2 == 1 
replace ecog = 3 if ecog3==1
replace ecog = 4 if ecog4==1


gen time_limit = 365
gen outcome = "progression"


gen insured = 0
replace insured = 1 if patientassistanceprogram == 1 | othergovernmentalinsurance == 1 | medicare == 1 | medicaid == 1 | commercialhealthplan ==1


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

gen histology = 2
replace histology = 0 if squamouscellcarcinoma == 1
replace histology = 1 if nonsquamouscellcarcinoma == 1


drop if timetocensor < 365 & (progressionoutcome == 0 | mortalityoutcome==0)

gen race = 5
replace race = 0 if white == 1
replace race = 1 if black == 1
replace race = 2 if asian == 1
replace race = 3 if hispanicrace == 1
replace race = 4 if otherrace == 1


gen gender = 2
replace gender = 0 if male == 1
replace gender = 1 if female == 1

keep if pdl1reported==1

gen smokinghistory = 2
replace smokinghistory = 0 if neversmoker == 1
replace smokinghistory = 1 if previoussmoker==1

gen uninsured = 0
replace uninsured = 1 if (noinsurance == 1 | selfpay==1) & medicare == 0 & medicaid == 0 & commercialhealthplan==0 & patientassistanceprogram== 0 & othergovernmentalinsurance== 0

gen otherinsurance = 0
replace otherinsurance = 1 if  (patientassistanceprogram== 1 | othergovernmentalinsurance== 1)


 replace progressionoutcome = 1 if progressiondays <= 365  & progressiondays>0
 replace mortalityoutcome = 1 if mortalitydays > 0 & mortalitydays <=365
 replace progressiondays = timetocensor if timetocensor < 365 & progressionoutcome==0
 replace mortalitydays = timetocensor if timetocensor < 365 & mortalityoutcome==0
 replace progressiondays = 365 if progressiondays > 365
 replace mortalitydays = 365 if mortalitydays > 365
  replace progressiondays = 365 if progressiondays == 0 | progressiondays > 365
  replace mortalitydays = 365 if mortalitydays == 0 | mortalitydays > 365
  
   drop if race == 3 | ecog==4
 teffects psmatch (progressionoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5,  nneighbor(1)
teffects psmatch (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5,  nneighbor(1)
teffects ipw (progressionoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5
teffects ipw (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5
teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1 >=0.5


 drop if pdl1>=0.5
 logit therapy_type $therapy_type_no_pdl1 
 predict yhat
 
 teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1) 
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1) 


graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density Pre-Matching") xtitle("Propensity Score") legend(label (1 "First-Line Chemotherapy") label(2 "First-Line IO Monotherapy"))
graph export "propensity_score_pre_matching_progression_reck.png", replace
graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1, ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Pre-Matching") ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "IO Monotherapy") ///
		label(2 "Chemoimmunotherapy"))
graph export "propensity_pre_match_hist_prog_reck.png", replace


pstest ${therapy_type_no_pdl1}, raw treated(therapy_type)

psmatch2 therapy_type  ${therapy_type_no_pdl1} , outcome(progressionoutcome) caliper(0.2) n(1) noreplacement 
count 

sum _weight
drop if _weight==.
drop if _support==0
count if _treated==1
count if _treated==0

pstest ${therapy_type_no_pdl1} , treated(therapy_type)

graph twoway (kdensity yhat if therapy_type==0) (kdensity yhat if therapy_type==1) ,ytitle("Propensity Score Density") xtitle("Propensity Score") legend(label (1 "Chemoimmunotherapy") label(2 "IO Monotherapy"))
graph export "propensity_post_match_reck.png", replace

graph twoway (histogram yhat if therapy_type==0, fcolor(blue%25) ///
		lcolor(none%0)) (histogram yhat if therapy_type==1,  ///
		fcolor(gray%25) lcolor(none%0)), xti("Propensity score") ///
		title("Distribution of Propensity Scores Post-Matching")  ///
		subtitle("by Therapy Type, Entire Dataset") legend(label(1 "Chemoimmunotherapy") ///
		label(2 "IO Monotherapy")) 
graph export "propensity_post_match_hist_prog_over_50.png", replace



 stset progressiondays, failure(progressionoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Progression-Free Survival") subtitle("PD-L1 < 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_over_50_pfs.png", replace
  sts test therapy_type, logrank
  


 stset mortalitydays, failure(mortalityoutcome)
 stci, by(therapy_type) rmean
 stcox therapy_type
 sts graph, by(therapy_type) title("Overall Survival") subtitle("PD-L1 < 50%") xtitle ("Survival Time from Treatment Initiation (Days)") ytitle("Proportion at Risk") legend(order(1 "Chemoimmunotherapy" 2 "IO Monotherapy")) 
 graph export "combo_io_post_match_over_50_os.png", replace
  sts test therapy_type, logrank
  
  
pstest ${therapy_type_no_pdl1}, treated(therapy_type)
