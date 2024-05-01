/*---------------------------------------------------------
HPR6070_Stata_Lab_4.do Instrumental Variables

	-> use ivregress and compare with OLS results
	-> use postestimation commands to assess assumptions
	-> bivariate probit
	-> exporting results in a table
	
Last Modified: March 27, 2024 (Elise Parrish)
	Thanks and credit to Kate Miller	
*---------------------------------------------------------*/

/*---------------------------------------------------------
		Setting up directories
-----------------------------------------------------------*/

clear
global path "/Users/vahluw/Documents/NSCLC_PDL1_Immunotherapy/"

cd "${path}"

set seed 6070
ssc install ivregress2 , replace
ssc install estout, replace
ssc install table1, replace
ssc install outreg2, replace

/*Import data directly from website using line of code below or
use the line of code: use "[insert file path]\StataLabIV.dta", clear*/
use http://fmwww.bc.edu/ec-p/data/wooldridge/card , clear

//If the above doesn't work-- the data are on canvas
	//use "card_data", clear

/*Question: what is the return to education on wages?

Card, D. "Using Geographic Variation in College Proximity to Estimate the Return to Schooling" https://davidcard.berkeley.edu/papers/geo_var_schooling.pdf

"This paper explores the use of college proximity as an exogenous determinant of school.

The marginal returns to education among children of less-education parents are as high and maybe higher than the rates of return estimated by conventional methods."

*/
/*---------------------------------------------------------
		Codebook
-----------------------------------------------------------*/
codebook

/*---------------------------------------------------------
		OLS
-----------------------------------------------------------*/

// Controlling for experience, experience^squared
reg lwage educ c.exper c.expersq i.black i.smsa i.south , robust 
	di exp(_b[educ])
	di (exp(_b[educ])-1)*100
estimates store ols

* why might our OLS estimate be biased?

/*---------------------------------------------------------
		Finding an instrument
-----------------------------------------------------------*/
//Distance to college (nearc4)
*conceptually, the argument is that proximity to a college has a direct effect on having more education but an indirect effect on wages

*what do key variables look like by proximity to nearest college?
table1 ,  vars ( ///
	lwage contn \ ///
	wage conts \ ///
	educ conts \ ///
	exper  conts \ ///
	black bin \ ///
	smsa bin \ ///
	south bin ) ///
	by (nearc4) format(%9.2f) cformat(%9.2f) test

//derfine covariates for later
global x "c.exper c.expersq i.black i.smsa i.south"

/*---------------------------------------------------------
		2SLS 
-----------------------------------------------------------*/

help ivregress
ivregress 2sls lwage (educ=nearc4) $x, first 
estimates store tsls

*----- Compare OLS and 2SLS -----*
esttab ols tsls , mti("OLS" "2SLS") se

	*recall also that 2SLS estimates the LATE

/*---------------------------------------------------------
		What is ivregress doing?
-----------------------------------------------------------*/

*----- If we manually attempted what the ivregress comand does-does it work? -----*
//first stage
reg educ i.nearc4 $x, robust
	predict educ_fitted //post estimation saving predicted values
//second stage
reg lwage c.educ_fitted $x, robust
estimates store manual

//Compare the coefficients and standard errors to ivregress output
esttab tsls manual , mti("ivregress" "by-hand") se order(educ educ_fitted)
	
/*---------------------------------------------------------
		Assessing assumptions 
-----------------------------------------------------------*/
help ivregress postestimation

*----- Instrument strength (or relevance) -----*
//Rule of thumb: Fstat>10
// Check the strength of the instrument --> is wage correlated with proximity to college
ivregress 2sls lwage  (educ=nearc4) $x, first robust
	estat firststage

reg educ i.nearc4 $x, robust //statistically significant coef in 1st stage
	test _b[1.nearc4] = 0

	
ivreghdfe lwage  (educ=nearc4) $x, first robust

	
/*
First-stage least squares regression: Linear regression of Instrument (proximity) --> Treatment (College degree)
Second-stage least squares regression: Linear regression of **PREDICTED** treatment (college degree) --> Outcome (wages)
*/	
	
*----- Instrument validity (or exogeneity and the exclusion restriction)-----*
//can only do a suggestive test here --mostly about thinking theoretically

//educ_fitted = predicted treatment based on instrument (predicted values)
//Put IV directly into the 2nd stage
/* HERE: We see that there is a significant relationship between our instrument (1.nearc4) and outcome (lwages), which means that maybe...our IV isn't perfect*/
	reg lwage c.educ_fitted i.nearc4 $x, robust



//No formal test for monotonicity either

*----- Other tests -----*	
ivregress 2sls lwage  (educ=nearc4) $x , first robust
	//Durbin-Wu-Hausman test (OLS vs 2SLS)
	*2SLS is less "efficient" than OLS so sometimes 
	* we want to see if 2SLS is necessary 
	* Null hypothesis: OLS estimator is consistent
	estat endogenous 
	
		//by hand:
		// Storing the residuals from the predicted regression
		reg educ i.nearc4 $x, robust
		predict v_hat, res
		
		//Then...determine relationship between outcome and 1) treatment and 2) residuals
		reg lwage educ v_hat $x, robust
		//Test the coefficient on the residuals
		//p-value = 0.2051
		test v_hat
			//The null hypothesis is the same as above

	
//can run this many many times using bootstrap
 cap program drop tsri_ex
program define tsri_ex, eclass
	reg educ i.nearc4 $x , robust
	predict res , res
	reg lwage c.educ c.res $x , robust
	drop res
end

bootstrap, reps(100): tsri_ex

	 
/*---------------------------------------------------------
		Bivariate probit (eg of a nonlinear case)
-----------------------------------------------------------*/

//consider a case with a binary treatment variable and binary outcome)
gen educ_bin = 0 if educ < 14
replace educ_bin = 1 if educ >= 14 & educ < .

gen wage_bin = 0 if wage < 540
replace wage_bin = 1 if wage >= 540 & wage < .


help biprobit
biprobit (wage_bin = educ_bin $x ) (educ_bin = i.nearc4 $x )


/*---------------------------------------------------------
		Exporting results from multiple models in one table
-----------------------------------------------------------*/
//Output Model Results
ivregress2 2sls lwage (educ=nearc4) $x, first 
est restore first
outreg2 using "Results.xls", replace excel dec(3)
est restore second
outreg2 using "Results.xls", append excel dec(3)

