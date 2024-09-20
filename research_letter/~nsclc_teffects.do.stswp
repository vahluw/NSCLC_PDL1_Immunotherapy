/*  Chemo therapy vs. first-line IO monotherapy Kaplan-Meier (PDL1 and non-PDL1) */



//global therapy_type_no_pdl1 "i.sdh braf kras i.ecog i.histology i.race   hispanicethnicity  diagnosisyear ageatdiagnosis i.gender  daysfromadvanceddiagnosistotreat   medicare selfpay medicaid commercialhealthplan otherinsurance uninsured  i.smokinghistory  academicmedicalcenter diabetes bonemetastases brainmetastases   antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen "

global therapy_type_no_pdl1 "i.sdh braf kras i.ecog i.histology i.race diagnosisyear ageatdiagnosis i.gender  daysfromadvanceddiagnosistotreat  i.smokinghistory  academicmedicalcenter antiinfectiveusepriortotreatment glucocorticoidusepriortotreatmen  medicare selfpay medicaid commercialhealthplan otherinsurance uninsured "


global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/research_letter/"
cd "${path}"
set scheme cleanplots

import delimited "all_data.csv", clear 

gen ecog = 5
replace ecog = 0 if ecog0==1
replace ecog = 1 if ecog1==1
replace ecog = 2 if ecog2==1 
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

drop if timetocensor < 182 & (progressionoutcome == 0 | mortalityoutcome==0)

replace progressionoutcome = 1 if progressiondays <= 182  & progressiondays>0
replace mortalityoutcome = 1 if mortalitydays <= 182  & mortalitydays>0
replace progressionoutcome = 0 if progressiondays > 182  | progressiondays==0
replace mortalityoutcome = 0 if mortalitydays > 182  | mortalitydays==0
replace progressiondays = 182 if progressiondays == 0 | progressiondays > 182
replace mortalitydays = 182 if mortalitydays == 0 | mortalitydays > 182


drop if race == 3 | ecog==4

teffects ipwra (progressionoutcome pdl1) (therapy_type)
teffects ipwra (mortalityoutcome pdl1) (therapy_type)
teffects ipwra (progressionoutcome pdl1) (therapy_type) if pdl1>=0.5
teffects ipwra (mortalityoutcome pdl1) (therapy_type) if pdl1>=0.5
teffects ipwra (progressionoutcome pdl1) (therapy_type) if pdl1<0.5 & pdl1>0
teffects ipwra (mortalityoutcome pdl1) (therapy_type) if pdl1<0.5 & pdl1>0
teffects ipwra (progressionoutcome pdl1) (therapy_type) if pdl1<0.5 
teffects ipwra (mortalityoutcome pdl1) (therapy_type) if pdl1<0.5

teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1)
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1)
teffects ipwra (progressionoutcome ) (therapy_type $therapy_type_no_pdl1)
teffects ipwra (mortalityoutcome) (therapy_type $therapy_type_no_pdl1)


teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1>=0.5
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1>=0.5
teffects ipwra (progressionoutcome ) (therapy_type $therapy_type_no_pdl1) if pdl1>=0.5
teffects ipwra (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1>=0.5

teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1<0.5
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1<0.5
teffects ipwra (progressionoutcome ) (therapy_type $therapy_type_no_pdl1) if pdl1<0.5
teffects ipwra (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1<0.5

teffects ipwra (progressionoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1<0.5 & pdl1>0
teffects ipwra (mortalityoutcome pdl1) (therapy_type $therapy_type_no_pdl1) if pdl1<0.5 & pdl1>0
teffects ipwra (progressionoutcome ) (therapy_type $therapy_type_no_pdl1) if pdl1<0.5 & pdl1>0
teffects ipwra (mortalityoutcome) (therapy_type $therapy_type_no_pdl1) if pdl1<0.5 & pdl1>0
