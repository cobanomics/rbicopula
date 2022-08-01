*! version 1.0.0 , 01aug2022
*! Author: Mustafa Coban, Institute for Employment Research (Germany)
*! Website: mustafacoban.de
*! Support: mustafa.coban@iab.de


*!***********************************************************!
*!     recursive bivariate copula estimation     		 	*!
*!***********************************************************!


/*	ESTIMATION	*/
	
program define rbicopula, eclass prop(ml_score svyb svyj svyr)
	version 11
	
	if replay(){
		if ("`e(cmd)'" != "rbicopula") error 301
		Display `0'
		exit `rc'
	}
	else{
			
		gettoken key: 0, parse(" =,[")
		
		*!	postestimation commands	
		if `"`key'"' == "margdec"{
					
			cap confirm variable `key'
			if !_rc{
				
				dis as err "Please rename the variable {bf:`key'} in the data."
				dis as err "This might lead to conflict with the postestimation command"  ///
							"{bf: rbicopula margdec}"
				exit 119
			}
			
			gettoken pfx 0 : 0, parse(" =,[")
			rbicopula_margdec `0'
		}
		else if `"`key'"' == "tmeffects"{
						
			cap confirm variable `key'
			if !_rc{
				
				dis as err 	"Please rename the variable {bf:`key'} in the data."
				dis as err 	"This might lead to conflict with the postestimation command"  ///
							"{bf: rbicopula tmeffects}"
				exit 119
			}
			
			gettoken pfx 0 : 0, parse(" =,[")
			rbicopula_tmeffects `0'
		}
		else{

			*!	rbiprobit estimation
			_vce_parserun rbicopula, mark(OFFset CLuster): `0'
			if "`s(exit)'" != ""{
				ereturn local cmdline "rbicopula `0'"
				exit
			}
			
			Estimate `0'
			ereturn local cmdline "rbicopula `0'"
		}
	}
end




program define Estimate, eclass
	
	
	*!	syntax based on biprobit and heckprobit
	gettoken dep1 0 : 0, parse(" =,[")		
	_fv_check_depvar `dep1'					
	
	tsunab dep1: `dep1'			
	rmTS `dep1'					
	
	confirm variable `r(rmTS)'	
	local dep1n 	`r(rmTS)'	
	local dep1eqn 	"`dep1'"
	
	*!	allow "=" after depvar
	gettoken equals rest : 0, parse(" =")
	if "`equals'" == "=" { 
		local 0 `"`rest'"' 
	}	
	
	
	
	*! parse indepvars
	syntax [varlist(numeric default=none ts fv)] [if] [in] [pw fw iw] , /*
			*/	ENDOGenous(string) [	/*	
			*/	COPula(string) noCONStant OFFset(varname numeric)	/*
			*/	Robust Cluster VCE(passthru) SCore(string) ITERate(passthru)	/*
			*/	moptobj(passthru) MLOpts(string) FROM(string) Level(cilevel)	/*
			*/	noSKIP LRMODEL noLOG	/*
			*/	* ]
	
	local fvops1 = "`s(fvops)'" == "true" | _caller() >= 11
	local tsops1 = "`s(tsops)'" == "true" | _caller() >= 11
			
	local indep1 	`varlist'
	local nocons1	`constant'
	local offset1	`offset'
	local option0 	`options'
	
	
	*!	parse depvar_en, indepvars_en and enopts
	EndogEq dep2 indep2 option2 : `"`endogenous'"'
	
	marksample 	touse
	markout		`touse' `dep1' `dep2' `indep1' `indep2' `offset1' `offset2'
	
	if inlist("`copula'", "", "gaussian"){
		macro drop	COPULAFUNC
		global 		COPULAFUNC = "gaussian"
		local copula "gaussian"
		local coptitle = strupper("`copula'")
	}
	else if inlist("`copula'","product","fgm","amh","plackett","clayton","frank","gumbel","joe"){
		macro drop	COPULAFUNC
		global 		COPULAFUNC = "`copula'"
		local coptitle = strupper("`copula'")
	}
	else{
		noi dis in red "option copula() incorrectly specified."
		error 198
	}

	
	if "`copula'" == "frank"{
		cap findfile integrate.ado
		if _rc == 601{
			noi dis as err "user-written command {bf:integrate} needed and not installed."
			noi dis as res "Type {stata ssc install integrate}"
			exit 601
		}
	}
	
	
	
	_get_diopts diopts option0, `option0'
	mlopts	stdopts, `option0'
	local	coll `s(collinear)'
	local	cns	`s(constraints)'
	
	
	if "`score'" != ""{
		local wdcnt: word count `score'
		if `wdcnt' == 1 & bsubstr("`score'",-1,1) == "*"{
			*!	User defined STUB*
			local score = bsubstr("`score'",1,length("`score'")-1)
			local score `score'1 `score'2 `score'3
			local wdcnt
		}
	}
	
	local scvar `"`score'"'
	
	if "`scvar'" != ""{
		local wdcnt: word count `scvar'
		if `wdcnt' != 3{
			noi dis in red "score() requires that 3 variables be specified"
			exit 198
		}
		confirm new var `scvar'
		local scvar1: word 1 of `scvar'
		local scvar2: word 2 of `scvar'
		local scvar3: word 3 of `scvar'
		tempvar	sc1 sc2 sc3
		local scopt "score(`sc1' `sc2' `sc3')"
	}
	

	local wtype `weight'
	local wtexp `"`exp'"'
	if "`weight'" != "" { 
		local wgt `"[`weight'`exp']"'  
	}
	
	if "`cluster'" != "" { 
		local clopt "cluster(`cluster')" 
	}
	
	_vce_parse, argopt(CLuster) opt(Robust oim opg) old: ///
				[`weight'`exp'], `vce' `clopt' `robust'	
	local cluster `r(cluster)'
	local robust `r(robust)'
	
	if "`cluster'" != ""{ 
		local clopt "cluster(`cluster')" 
	}
	
	if `"`robust'"' != "" {
		local crtype crittype("log pseudolikelihood")
	}

		
	if _caller() < 15{
		local parm delta:_cons
	}
	else{
		local parm /:delta
	}

	
	if inlist("`copula'","product"){
		local diparm_tet ""
	}
	else if inlist("`copula'","gaussian","fgm","amh"){
		local diparm_tet "diparm(delta, tanh label(theta))"
	}
	else if inlist("`copula'","plackett","clayton"){
		local diparm_tet "diparm(delta, exp label(theta))"
	}	
	else if inlist("`copula'","frank"){
		local diparm_tet "diparm(delta, function(@) derivative(1) label(theta))"
	}
	else if inlist("`copula'","gumbel","joe"){
		local diparm_tet "diparm(delta, function(1+exp(@)) derivative(exp(@)) label(theta))"
	}		
	     
	
	if "`copula'" == "product"{
		local diparm_tau ""
	}
	else if "`copula'" == "gaussian"{
		local diparm_tau	diparm(delta, function(2*asin(tanh(@))/_pi) ///
							derivative(2*sqrt(1-tanh(@)^2)/_pi) label("tau"))
	}
	else if "`copula'" == "fgm"{
		local diparm_tau	diparm(delta, function(2*tanh(@)/9) ///
							derivative(2*(1-tanh(@)^2)/9) label("tau"))
	}
	else if "`copula'" == "plackett"{
		local diparm_tau ""
	}
	else if "`copula'" == "clayton"{
		local diparm_tau	diparm(delta, function(exp(@)/(exp(@)+2)) ///
							derivative(2*exp(@)/((exp(@)+2)^2)) label("tau"))
	}
	else if "`copula'" == "frank"{
		local diparm_tau ""
	}
	else if "`copula'" == "gumbel"{
		local diparm_tau	diparm(delta, function(1-1/(1+exp(@))) ///
							derivative((1+exp(@))^(-2)*(exp(@))) label("tau"))
	}
	else if "`copula'" == "joe"{
		local diparm_tau diparm(delta, function(1 + (2/(2-(1+exp(@)))) * (digamma(2) - ///
				digamma(1 + (2/(1+exp(@)))))) ///
				derivative( ((2*digamma(2))/((2-(1+exp(@)))^2)) - ///
				((2*digamma(1 + (2/(1+exp(@))))) / ((2-(1+exp(@)))^2)) + ///
				((4*trigamma(1 + (2/(1+exp(@))))) / ((2-(1+exp(@)))*(1+exp(@))^2))) ///
				label("tau"))
	}
	else if "`copula'" == "amh"{
		local diparm_tau	diparm(delta, function((3*tanh(@)-2)/(3*tanh(@)) - ///
				(2/3)*(1-1/tanh(@))^2 *ln(1-tanh(@))) ///
				derivative(-(2/3)*(-2+tanh(@)+2*(1-1/tanh(@))*ln(1-tanh(@)))/(tanh(@)^2)*(1-tanh(@)^2)) ///
				label("tau"))
	}

	
	if "`log'" == ""{
		local log "noisily"
	}
	else{
		local log "quietly"
	}
	

	if "`level'" != ""{
		local level "level(`level')"
	}
	

	_rmcoll i.`dep2' `indep1' `wgt' if `touse' , `nocons1' `coll'
	local indep1 "`r(varlist)'"
	
	_rmcoll `indep2' `wgt' if `touse' , `nocons2' `coll'
	local indep2 "`r(varlist)'"
	
	
	if "`offset1'" != "" { 
		local offo1 "offset(`offset1')"
	}
	if "`offset2'" != "" { 
		local offo2 "offset(`offset2')" 
	}
	
	
	if "`lrmodel'" != ""{
		_check_lrmodel , `skip' `constant' constraints(`cns')	/*
					*/	options(`clopt' `robust' `nocons1' `nocons2')	/*
					*/	indep(`indep1' `indep2')
		local skip	noskip
	}
	local skip = cond("`skip'" != "", "", "skip")
	
	
	
	*!	Taken from "bicop" and modified a bit	
	qui{
		count	if `touse'
		local N = r(N)
		if `N' == 0{
			error 2000		
		}

		tab `dep1n'	if `touse'
		local nthr1 = r(r)

		tab `dep2n'	if `touse'
		local nthr2 = r(r)
		
		if `nthr1' != 2{
			if `nthr1' == 1{
				dis in red "There is no variation in `dep1'"
				exit 2000
			}
			else if `nthr1' > 2{
				dis in red "There are more than two groups in `dep1'"
				exit 2000			
			}
		}
		
		if `nthr2' != 2{
			if `nthr2' == 1{
				dis in red "There is no variation in `dep2'"
				exit 2000
			}
			else if `nthr2' > 2{
				dis in red "There are more than two groups in `dep2'"
				exit 2000				
			}
		}
	}
	
	
	*!	0/1 values for depvar and depvar_en	
	qui: levelsof `dep1n'
	if "`r(levels)'" != "0 1"{
		dis in green "{bf:`dep1'} does not vary; remember:"
        dis in green "0 = negative outcome, 1 = positive outcome"
		exit 2000
	}
	
	qui: levelsof `dep2n'
	if "`r(levels)'" != "0 1"{
		dis in green "{bf:`dep2'} does not vary; remember:"
        dis in green "0 = negative outcome, 1 = positive outcome"
		exit 2000
	}
	
	
	if "`from'" == ""{

		`log' dis in green _n "Univariate Probits for starting values"
	
		*!	eq.1: univariate probit (taken from biprobit)
		`log' dis in green _n "Fitting comparison outcome equation:"
		cap `log' probit `dep1' `indep1' `wgt' if `touse' , /*
					*/	`offo1' `nocons1' nocoef asis `stdopts'	`crtype' /*
					*/	iter(`=min(1000,c(maxiter))') 
				
		if _rc == 0{
			tempname cb1
			mat `cb1' 	= e(b)
			local ll_1	= e(ll)
			mat coleq `cb1' = `dep1eqn'
		}
	

		*!	eq.2: univariate probit (taken from biprobit)
		`log' dis in green _n "Fitting comparison treatment equation:"	
		cap `log' probit `dep2' `indep2' `wgt' if `touse' , /*
					*/	`offo2' `nocons2' nocoef asis `stdopts'	`crtype' /*
					*/	iter(`=min(1000,c(maxiter))')

		local ll_str = e(crittype)
											
		if _rc == 0{
			tempname cb2
			mat `cb2'	= e(b)
			local ll_2 	= e(ll)
			mat coleq `cb2' = `dep2eqn'
		}
	
	
		if "`ll_1'" == ""{
			local ll_1 = 0
		}
		if "`ll_2'" == ""{
			local ll_2 = 0
		}
		
		local ll_p = `ll_1' + `ll_2'
	
		*!	stack coefficient estimates
		if "`cb1'`cb2'" != ""{
			tempname from
			
			if "`cb1'" != "" & "`cb2'" != ""{
				mat `from' = `cb1', `cb2'
				`log' dis in green _n "Comparison:    `ll_str' = " in yellow %10.0g `ll_p'
			}
			else if "`cb1'" != ""{
				mat `from' = `cb1'
				local lrtest "nolrtest"
			}
			else{
				mat `from' = `cb2'
				local lrtest "nolrtest"
			}
			local getvals 1
		}
	}
	
	
	if "`nocons1'`nocons2'" != ""{
		local skip "skip"
	}
	if "`indep1'`indep2" == ""{
		local skip = "skip"
	}
	if "`robust'" != ""{
		local skip = "skip"
	}
	
	
	*!	constants-only model
	tempname a0
	if "`skip'" == ""{
		
		*!	eq.1: univariate probit (taken from biprobit)
		cap qui probit `dep1' i.`dep2' `wgt' if `touse' , /*
					*/	`offo1' asis `stdopts' iter(`=min(1000,c(maxiter))') 
				
		if _rc == 0{
			tempname ccb1
			mat `ccb1' 	= e(b)
			mat coleq `ccb1' = `dep1eqn'
		}		
	
		*!	eq.2: univariate probit (taken from biprobit)
		cap qui probit `dep2' `wgt' if `touse' , /*
					*/	`offo1' asis `stdopts' iter(`=min(1000,c(maxiter))') 
				
		if _rc == 0{
			tempname ccb2
			mat `ccb2' 	= e(b)
			mat coleq `ccb2' = `dep2eqn'
		}	
		
		*!	stack coefficient estimates
		tempname a
		if "`copula'" == "product"{
			mat `a' = (0)
		}
		else{
			mat `a' = (-.3)		//	biprobit uses (-.3)
		}
				
		mat colnames `a' = `parm'
		
		if "`cb1'`cb2'" != ""{
			tempname from0
			
			if "`cb1'" != "" & "`cb2'" != ""{
				mat `from0' = `ccb1', `ccb2'
			}
			else if "`ccb1'" != ""{
				mat `from0' = `ccb1'
			}
			else{
				mat `from0' = `ccb2'
			}
		}
		
		mat `from0' = `from0',`a'
		
		
		
		*!	constans-only model estimation
		`log' dis in green _n "Fitting constant-only model"	
		
		#d ;
		`log' ml model lf1 rbicopula_lf1
			(`dep1eqn': `dep1' = i.`dep2', `nocons1' `offo1')
			(`dep2eqn': `dep2' = , `nocons2' `offo2')
			/delta
			if `touse' `wgt'
			,
			maximize init(`from0') search(off) wald(0) `iterate' `mlopts'	
			`stdopts' nooutput missing nopreserve collinear	`level'	
			`diparm_tet' `diparm_tau' `lrtest' `diopts'	nocnsnotes	`crtype' negh
			;
		#d cr
	
	
		*!	initial value for delta (taken from biprobit)	
		local cont "continue"
		if "`getvals'" != ""{
			mat `a0' = e(b)
			mat `a0' = `a0'[1,5]
			mat colnames `a0' = `parm'
			mat `from' = `from',`a0'
		}
	}
	else{
		if "`getvals'" != ""{
			if "`copula'" == "product"{
				mat `a0' = (0)
			}
			else{
				mat `a0' = (0.01)		// heckmancopula uses (0.01)
			}
			
			mat colnames `a0' = `parm'
			mat `from' = `from', `a0'
		}	
	}
	

	if "`cns'" != "" | "`skip'" != ""{
		local cont `cont' wald(2)
	}
	
	

	*!	full model estimation
	`log' dis in green _n "Fitting full model:"	
	local title "Recursive Bivariate Copula Regression (Copula: `coptitle')"
	
	#d ;
	`log' ml model lf1 rbicopula_lf1
			(`dep1eqn': `dep1' = `indep1', `nocons1' `offo1')
			(`dep2eqn': `dep2' = `indep2', `nocons2' `offo2')
			/delta
			if `touse' `wgt'
			,
			maximize init(`from') search(off) `cont' `robust' `clopt' `scopt'			
			`iterate' `moptobj' `mlopts' `stdopts' nooutput	missing	nopreserve collinear		
			title(`title') `level' `diparm_tet' `diparm_tau' `lrtest' `diopts' negh
			;
	#d cr
	
	
	
	*!	constraints matrix
	tempname junk
	capture matrix `junk' = e(Cns)
	if !_rc{
		local hascns hascns	
	}
	
	if "`scvar'" != ""{
		rename `sc1' `scvar1'
		rename `sc2' `scvar2'
		rename `sc3' `scvar3'
		ereturn local scorevars `scvar'
	}
	

	*!	estimate of dependence parameter	
	if inlist("`copula'","product"){
		ereturn scalar theta = 0
	}
	else if inlist("`copula'","gaussian","fgm","amh"){
		ereturn scalar theta = tanh(_b[/delta])
	}
	else if inlist("`copula'","plackett","clayton"){
		ereturn scalar theta = exp(_b[/delta])
	}	
	else if inlist("`copula'","frank"){
		ereturn scalar theta = _b[/delta]
	}
	else if inlist("`copula'","gumbel","joe"){
		ereturn scalar theta = 1 + exp(_b[/delta])
	}	
	
	
	*!	estimate of kendall's tau
	if inlist("`copula'","product"){
		ereturn scalar tau = 0
	}
	else if "`copula'" == "gaussian"{
		ereturn scalar tau = 2*asin(`e(theta)')/_pi
	}	
	else if "`copula'" == "fgm"{
		ereturn scalar tau = 2*`e(theta)'/9
	}	
	else if "`copula'" == "clayton"{
		ereturn scalar tau = `e(theta)'/(`e(theta)'+2)
	}	
	else if "`copula'" == "frank"{
		mata: Tau_Value(`e(theta)')
		ereturn scalar tau = tau
	}
	else if "`copula'" == "gumbel"{
		ereturn scalar tau = (`e(theta)'-1)/`e(theta)'
	}	
	else if "`copula'" == "joe"{
		ereturn scalar tau = (1 + (2/(2-`e(theta)')) * (digamma(2) - digamma(1 + (2/`e(theta)'))))
	}	
	else if "`copula'" == "amh"{
		ereturn scalar tau = ((3*`e(theta)'-2)/(3*`e(theta)') - ///
				(2/3)*(1-1/`e(theta)')^2 *ln(1-`e(theta)'))
	}	

	
	
	if "`copula'" != "product"{
		if "`ll_p'" != "" & "`robust'`lrtest'`hascns'" == ""{
			*!	likelihood-ratio test
			ereturn scalar	ll_c 	= `ll_p'
			ereturn scalar	chi2_c	= abs(-2*(e(ll_c) - e(ll)))	//	
			ereturn local	chi2_ct "LR"
		}
		else{
			*!	wald test			
			if inlist("`copula'","gaussian","fgm","amh"){
				local test_tet (tanh(_b[/delta]) = 0)
			}
			else if inlist("`copula'","plackett"){
				local test_tet (exp(_b[/delta]) = 1)
			}	
			else if inlist("`copula'","clayton"){
				local test_tet (exp(_b[/delta]) = 0)
			}	
			else if inlist("`copula'","frank"){
				local test_tet (_b[/delta] = 0)
			}
			else if inlist("`copula'","gumbel","joe"){
				local test_tet ((1 + exp(_b[/delta])) = 1)
			}			
			
			cap qui testnl `test_tet'		
			ereturn scalar chi2_c = r(chi2)	
			ereturn local chi2_ct "Wald"
		}
	}


	*!	stored results
	ereturn scalar k_aux = 1	
	ereturn scalar k_eq_model = 2		
	ereturn local marginsok		"default P11 P10 P01 P00 PMARG1 PMARG2 PCOND1 PCOND2 PCOND10 XB1 XB2 PMARGCOND1"
	ereturn local marginsnotok	"STDP1 STDP2"
	ereturn hidden local marginsderiv 	"default P11 P10 P01 P00 PMARG1 PMARG2 PCOND1 PCOND2 PCOND10 PMARGCOND1"
	ereturn local option2	"`option2'"
	ereturn local option0	"`option0'"
	ereturn local copula	"`copula'"
	ereturn local predict	"rbicopula_p"
	ereturn local cmd 		"rbicopula"
	
	*!	display results
	Display, `level' `diopts'
	exit `e(rc)'
end




/* Debye function: integration */

mata
	real rowvector debye(real rowvector x)
	{
		return(x:/expm1(x))
	}
end


/* tau formula for frank and joe copula */
mata 
void Tau_Value(
	real scalar tet
)
{
	integ = integrate(&debye(),0,tet)
	tau = 1 - (4/tet)*(1 - (1/tet)*integ)
	st_numscalar("tau",tau)
}

end


/* DISPLAY */

program define Display
	
	version 11
	syntax [, Level(cilevel) *]
	_get_diopts diopts, `options'
	
	
	if "`e(copula)'" == "product"{ 
		ml display, level(`level') nofootnote neq(2) `diopts'
	}
	else if inlist("`e(copula)'","frank"){
		ml display, level(`level') nofootnote `diopts' diparm(__lab__, label("tau") value(`e(tau)'))
	}
	else{
		ml display, level(`level') nofootnote `diopts'
	}

	DispLr
	_prefix_footnote

end


*!	taken from biprobit
program define DispLr
	
	version 11
	
	if "`e(ll_c)'`e(chi2_c)'" == "" {
		exit
	}
	
	local chi : di %8.0g e(chi2_c)
	local chi = trim("`chi'")
	
	if inlist("`e(copula)'","plackett","gumbel","joe"){
		local indepvalue 1
	}
	else{
		local indepvalue 0
	}
	
	
	if "`e(ll_c)'"=="" {
		if inlist("`e(copula)'","clayton","gumbel","joe"){
			local chi2pvalue (1-(.5*1+.5*chi2(1,e(chi2_c))))	// taken from heckmancopula
		}
		else{
			local chi2pvalue (chiprob(1,e(chi2_c)))
		}
		di in green "Wald test of theta=`indepvalue': " ///
			in green "chi2(" in ye "1" in gr ") = " ///
			in ye "`chi'" ///
			in green _col(59) "Prob > chi2 = " in ye %6.4f ///
			`chi2pvalue'
		exit
	}
	
	local chi2pvalue chiprob(1,e(chi2_c))
	di in green "LR test of theta=`indepvalue': " ///
		in green "chi2(" in ye "1" in gr ") = " in ye `chi' ///
		in green _col(59) "Prob > chi2 = " in ye %6.4f ///
		`chi2pvalue'
end




	
	
/* AUXILIARY SUB-PROGRAMS */	

*!	taken from biprobit
program define rmTS, rclass
	
	version 11
	
	local tsnm = cond( match("`0'", "*.*"),  		/*
			*/ bsubstr("`0'", 			/*
			*/	  (index("`0'",".")+1),.),     	/*
			*/ "`0'")

	return local rmTS `tsnm'
end




*!	parse endog() option
program define EndogEq
	
	version 11
	
	args dep2 indep2 option2 colon endog_eq	
	
	gettoken dep rest: endog_eq, parse(" =")		
	_fv_check_depvar `dep'						
		
	tsunab dep: `dep'											
	rmTS `dep'				
	confirm variable `r(rmTS)'	
	
	c_local dep2n	`r(rmTS)'	
	c_local `dep2' 	`dep'		
	c_local dep2eqn `dep'		
	
	
	*!	allow "=" after depvar_en
	gettoken equals 0 : rest, parse(" =")
	if "`equals'" != "=" { 
		local 0 `"`rest'"' 	
	}	
	
		
	
	*!	parse indepvar_en (based on biprobit and heckprobit) 	
	syntax [varlist(numeric default=none ts fv)] [, /*
			*/	noCONStant OFFset(varname numeric) ]
	
	
	local fvops2 = "`s(fvops)'" == "true" | _caller() >= 11
	local tsops2 = "`s(tsops)'" == "true" | _caller() >= 11
			
	
	c_local `indep2'	`varlist'
	c_local offset2		`offset'
	c_local nocons2		`constant'
end
