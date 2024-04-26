global covariates "i.race i.ecog i.histology i.stage  i.practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan"

import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/all_data.csv", clear
gen therapy_type = -1
replace therapy_type = 0 if first_line_chemo == 1
replace therapy_type = 1 if io_mono == 1
gen first_line = 0
replace first_line = 1 if therapy_type ==1
rdplot first_line pdl1, c(0.5)
binscatter first_line pdl1, rd(0.5) yti("Probability of IO Monotherapy Treatment") xti("Recentered PD-L1") 
	
graph export "discontinuity_treat.png", replace 

graph twoway (hist pdl1) , xline(0.5, lcolor(red))

kdensity pdl1 , xline(0.5)

//Plotting all, testing only within the optimal bandwidth estimated
rddensity pdl1 , pl c(0.5)

/* Actually do statistical analysis for RD without nivolumab for IO vs chemo */

import delimited "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/all_data.csv", clear
gen therapy_type = -1
replace therapy_type = 0 if first_line_chemo == 1
replace therapy_type = 1 if io_mono == 1
keep if pdl1_given==1 & therapy_type>=0
gen first_line = 0
replace first_line = 1 if therapy_type ==1
gen over_threshold = 0
replace over_threshold = 1 if pdl1>=0.5
//rdwinselect pdl1 race ecog histology stage  practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan, cutoff(0.5) nwindows(20)
rdrandinf progression_12 pdl1, cutoff(0.5) fuzzy(first_line) kernel(triangular) covariates(ecog  practice_type diag_year age_at_diagnosis other_no_insurance  self_pay pt_assistance other_gov_insurance medicare medicaid commercial_health_plan)
rddensity pdl1, c(0.5) // check sorting/bunching assumption
