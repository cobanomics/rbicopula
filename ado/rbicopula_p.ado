*! version 1.1.0 , 10oct2022
*! Author: Mustafa Coban, Institute for Employment Research (Germany)
*! Website: mustafacoban.de
*! Support: mustafa.coban@iab.de


/****************************************************************/
/*    			 rbicopula prediction							*/
/****************************************************************/


program define rbicopula_p, eclass
	
	version 11
	syntax [anything] [if] [in] [, SCores *]
	
	if ("`e(cmd)'" != "rbicopula"){
		error 301
		dis in red "rbicopula was not the last command"
	}
	
	if `"`scores'"' != ""{
		ml score `0'	
		exit
	}
	
	local myopts P11 P10 P01 P00 PMARG1 PMARG2 XB1 XB2
	local myopts `myopts' PCOND1 PCOND2 PCOND10 PMARGCOND1 STDP1 STDP2
	local myopts `myopts' d1(string) d2(string)			
	
	_pred_se "`myopts'" `0'
	
	if (`s(done)') exit			
	local vtyp 	`s(typ)'			
	local varn 	`s(varn)'			
	local 0	`"`s(rest)'"'				
	
	
	*!	parse predict
	syntax [if] [in] [, `myopts' noOFFset]
	
	local type `p11'`p10'`p01'`p00'`pmarg1'`pmarg2'`xb1'`xb2'
	local type `type' `pcond1'`pcond2'`pcond10'`pmargcond1'`stdp1'`stdp2'
	
	tokenize `e(depvar)'
	local dep1 `1'
	local dep2 `2'
	
	tsunab dep1: `dep1'
	tsunab dep2: `dep2'
				
	rmTS `dep1'
	confirm variable `r(rmTS)'
	local dep1n 	`r(rmTS)'	
	rmTS `dep2'
	confirm variable `r(rmTS)'
	local dep2n 	`r(rmTS)'
	
	local copfc = "`e(copula)'"

	
	*!	linear index xb1
	if "`type'" == "xb1"{
		local pred "Linear Prediction of `dep1'"
		
		_predict `vtyp' `varn' `if' `in', eq(#1) `offset'
		label var `varn' "`pred'"
		
		exit
	}	
	
	
	*!	linear index xb2
	if "`type'" == "xb2"{
		local pred "Linear Prediction of `dep2'"
		
		_predict `vtyp' `varn' `if' `in', eq(#2) `offset'
		label var `varn' "`pred'"
		
		exit
	}	
	
	
	*!	standard error of linear index xb1
	if "`type'" == "stdp1"{
		local pred "S.E. of Linear Prediction of `dep1'"
		
		_predict `vtyp' `varn' `if' `in', stdp eq(#1) `offset'
		label var `varn' "`pred'"
		
		exit
	}	
	
	
	*!	standard error of linear index xb2
	if "`type'" == "stdp2"{
		local pred "S.E. of Linear Prediction of `dep2'"
		
		_predict `vtyp' `varn' `if' `in', stdp eq(#2) `offset'
		label var `varn' "`pred'"
		
		exit
	}	
	
	
	
	*!	dependence parameter
	tempname del
	
	if `:colnfreeparms e(b)' {
		scalar `del' = _b[/delta]			//	version 15 upwards
	}
	else {
		scalar `del' = [delta]_b[_cons]		//	older versions
	}	
	
	
	*!	first and second derivatives of dependence parameter
	tempname	tet d1tetd1del d2tetd2del
	
	if inlist("`copfc'","product"){
		scalar `tet' 		= 0
		scalar `d1tetd1del' = 0
		scalar `d2tetd2del' = 0
	}
	else if inlist("`copfc'","gaussian","fgm","amh"){
		scalar `tet' 		= (exp(2*`del')-1) / (1+exp(2*`del'))
		scalar `d1tetd1del' = (4 * exp(2*`del')) / ( (1+exp(2*`del'))*(1+exp(2*`del')) )
		scalar `d2tetd2del' = 8*exp(2*`del') * (1-exp(2*`del')) ///
							/ ( (1+exp(2*`del'))*(1+exp(2*`del'))*(1+exp(2*`del')) )
	}
	else if inlist("`copfc'","plackett","clayton"){
		scalar `tet' 		= exp(`del')
		scalar `d1tetd1del' = exp(`del')
		scalar `d2tetd2del' = exp(`del')
	}	
	else if inlist("`copfc'","frank"){
		scalar `tet' 		= `del'
		scalar `d1tetd1del' = 1
		scalar `d2tetd2del' = 0
	}
	else if inlist("`copfc'","gumbel","joe"){
		scalar `tet'		= (1 + exp(`del'))
		scalar `d1tetd1del' = exp(`del')
		scalar `d2tetd2del' = exp(`del')
	}	
	
	
	

	if !inlist("`type'","pmarg1","pmarg2","xb","zg","stdp1","stdp2"){
		*!	linear index for probabilities
		tempvar xb zg dep2orig
		
		qui{
			clonevar `dep2orig' = `dep2n'
		
			if inlist("`type'","","p11","p01","pcond1", "pcond2"){
				replace `dep2n' = 1
			}
			else if inlist("`type'","p10","p00","pcond10"){
				replace `dep2n' = 0
			}
			
			_predict double `xb' `if' `in', eq(#1) `offset'
			_predict double `zg' `if' `in', eq(#2) `offset'
			
			replace `dep2n' = `dep2orig'
		}
	
	
		*!	algebraic signs
		tempname s t q1 q2
		
		if inlist("`type'","p11","pcond1","pcond2","pmargcond1"){
			scalar `s'	= 1
			scalar `t'	= 1
		}
		else if inlist("`type'","p10","pcond10"){
			scalar `s' = 1
			scalar `t' = 0
		}
		else if "`type'" == "p01"{
			scalar `s' = 0
			scalar `t' = 1
		}
		else if "`type'" == "p00"{
			scalar `s' = 0
			scalar `t' = 0	
		}
		
		scalar `q1' = 2*`s' - 1
		scalar `q2' = 2*`t' - 1

		
		
		*!	shortcuts of marginals
		local u			(normal(-`xb'))
		local v			(normal(-`zg'))
		local d1ud1xb	(-normalden(-`xb'))
		local d1vd1zg	(-normalden(-`zg'))
		local d2ud2xb	(`xb'*normalden(-`xb'))
		local d2vd2zg	(`zg'*normalden(-`zg'))
	}
	

	
	*!	predictions of copulas
	if `"`d1'`d2'"' == ""{
		
		if "`copfc'" == "product"{
			local cop	(`u'*`v')
		}
		if "`copfc'" == "gaussian"{
			local cop 	(binormal(-`xb', -`zg', `tet'))
		}
		if "`copfc'" == "fgm"{
			local util	(1 - `u')
			local vtil	(1 - `v')
			local cop	(`u'*`v'*(1 + `tet'*`util'*`vtil'))
		}
		if "`copfc'" == "plackett"{
			local tettil	(`tet' - 1)
			local r			(1 + `tettil'*(`u'+`v'))
			local A			(`r'*`r' - 4*`u'*`v'*`tet'*`tettil')
			local cop		((`r' - sqrt(`A')) / (2*`tettil'))
		}
		if "`copfc'" == "clayton"{
			local tettil	(1/`tet')
			local A			(`u'^(-`tet') + `v'^(-`tet') - 1)
			local cop		(`A'^(-`tettil'))
		}
		if "`copfc'" == "frank"{
			local tettil	(1/`tet')
			local util		(exp(-`tet'*`u'))
			local vtil		(exp(-`tet'*`v'))
			local eta		(exp(-`tet'))
			local A			(1 + ((`util'-1)*(`vtil'-1) / (`eta'-1)))
			local cop		(-`tettil'*ln(`A'))
		}
		if "`copfc'" == "gumbel"{
			local tettil	(1/`tet')
			local util		(-ln(`u'))
			local vtil		(-ln(`v'))
			local A			(`util'^(`tet') + `vtil'^(`tet'))
			local cop		(exp(-(`A'^(`tettil'))))
		}
		if "`copfc'" == "joe"{
			local tettil	(1/`tet')
			local util		(1 - `u')
			local vtil		(1 - `v')
			local A			((`util'^(`tet') + `vtil'^(`tet') - (`util'^(`tet'))*(`vtil'^(`tet'))))
			local cop		(1 - `A'^(`tettil'))
		}
		if "`copfc'" == "amh"{
			local util		(1 - `u')
			local vtil		(1 - `v')
			local A			(1 - `tet'*`util'*`vtil')
			local cop 		(`u'*`v')/`A'
		}	
	}
	
	
	
	*!	predictions of copulas for margins
	if `"`d1'`d2'"' != ""{
		
		if "`copfc'" == "product"{
			local cop		(`u'*`v')
		}
		if "`copfc'" == "gaussian"{
			local cop 		(binormal(-`xb', -`zg', `tet'))
		}
		if "`copfc'" == "fgm"{
			local util		(1 - `u')
			local vtil		(1 - `v')
			local cop		(`u'*`v'*(1 + `tet'*`util'*`vtil'))
		}
		if "`copfc'" == "plackett"{
			local tettil	(`tet' - 1)
			local r			(1 + `tettil'*(`u'+`v'))
			local A			(`r'*`r' - 4*`u'*`v'*`tet'*`tettil')
			local cop		((`r' - sqrt(`A')) / (2*`tettil'))
		}
		if "`copfc'" == "clayton"{
			local tettil	(1/`tet')
			local A			(`u'^(-`tet') + `v'^(-`tet') - 1)
			local cop		(`A'^(-`tettil'))
		}
		if "`copfc'" == "frank"{
			local tettil	(1/`tet')
			local util		(exp(-`tet'*`u'))
			local vtil		(exp(-`tet'*`v'))
			local eta		(exp(-`tet'))
			local A			(1 + ((`util'-1)*(`vtil'-1) / (`eta'-1)))
			local cop		(-`tettil'*ln(`A'))
		}
		if "`copfc'" == "gumbel"{
			local tettil	(1/`tet')
			local util		(-ln(`u'))
			local vtil		(-ln(`v'))
			local A			(`util'^(`tet') + `vtil'^(`tet'))
			local cop		(exp(-(`A'^(`tettil'))))
		}
		if "`copfc'" == "joe"{
			local tettil	(1/`tet')
			local util		(1 - `u')
			local vtil		(1 - `v')
			local A			((`util'^(`tet') + `vtil'^(`tet') - (`util'^(`tet'))*(`vtil'^(`tet'))))
			local cop		(1 - `A'^(`tettil')
		}
		if "`copfc'" == "amh"{
			local util		(1 - `u')
			local vtil		(1 - `v')
			local A			(1 - `tet'*`util'*`vtil')
			local cop 		(`u'*`v')/`A'
		}	
		
		
		*!	shortcuts
		if "`copfc'" == "gaussian"{
			
		}
		if "`copfc'" == "plackett"{
			local d1Ad1tet	(2*`r'*(`u'+`v') - 4*`u'*`v'*(`tettil'+`tet'))
		}
		if "`copfc'" == "clayton"{
			local d1Ad1tet	(-`u'^(-`tet')*ln(`u') - `v'^(-`tet')*ln(`v'))
			local d2Ad2tet	(`u'^(-`tet')*(ln(`u'))^2 + `v'^(-`tet')*(ln(`v'))^2)
		}
		if "`copfc'" == "frank"{
			local d1Ad1tet	((1/(`eta'-1)) * (`util'*`u'*(1-`vtil') + `vtil'*`v'*(1-`util') ///
							+ `eta'*(`A'-1))) 
			local d2Ad2tet	((1/(`eta'-1)) * (`util'*`vtil'*((`u'+`v')^2) - (`u'^2*`util') - ///
							(`v'^2*`vtil') -`eta'*(`A'- (2*`d1Ad1tet') - 1)))
		}
		if "`copfc'" == "gumbel"{
			local d1Ad1tet	(`util'^(`tet')*ln(`util') + `vtil'^(`tet')*ln(`vtil'))
			local d2Ad2tet	(`util'^(`tet')*ln(`util')^2 + `vtil'^(`tet')*ln(`vtil')^2)
		}
		if "`copfc'" == "joe"{
			local d1Ad1tet	(`util'^(`tet')*ln(`util')*(1-`vtil'^(`tet')) ///
							+ `vtil'^(`tet')*ln(`vtil')*(1-`util'^(`tet')))
			local d2Ad2tet 	(ln(`util')*(`d1Ad1tet' - `vtil'^(`tet')*ln(`vtil')) ///
							+ ln(`vtil')*(`d1Ad1tet' - `util'^(`tet')*ln(`util')))
		}
		
		
		
		*!	first derivatives of copulas
		if "`copfc'" == "product"{
			local d1copd1u 		(`v')
			local d1copd1v		(`u')
			local d1copd1tet	(0)
		}
		if "`copfc'" == "gaussian"{
			
		}
		if "`copfc'" == "fgm"{
			local d1copd1u 		(`v' + `tet'*`vtil'*`v' - 2*`u'*`tet'*`vtil'*`v')
			local d1copd1v		(`u' + `tet'*`util'*`u' - 2*`v'*`tet'*`util'*`u')
			local d1copd1tet	(`util'*`vtil'*`u'*`v')		
		}
		if "`copfc'" == "plackett"{
			local d1copd1u 		(.5 + ((2*`v'*`tet' - `r')/(2*sqrt(`A'))))
			local d1copd1v		(.5 + ((2*`u'*`tet' - `r')/(2*sqrt(`A'))))
			local d1copd1tet	((1/(`tettil'*sqrt(`A'))) * (`cop' - `u'*`v'))	
		}
		if "`copfc'" == "clayton"{
			local d1copd1u 		(`A'^(-`tettil'-1) * `u'^(-`tet'-1))
			local d1copd1v		(`A'^(-`tettil'-1) * `v'^(-`tet'-1))
			local d1copd1tet	( `tettil'*`cop'*(`tettil'*ln(`A') - (`d1Ad1tet'/`A')))	
		}
		if "`copfc'" == "frank"{
			local d1copd1u 		((`util'/`A') * (`vtil'-1)/(`eta'-1))
			local d1copd1v		((`vtil'/`A') * (`util'-1)/(`eta'-1))
			local d1copd1tet	((`tettil'^2)*ln(`A') - (`tettil'/`A')*`d1Ad1tet')	
		}
		if "`copfc'" == "gumbel"{
			local d1copd1u 		(`cop'*(`A'^(`tettil'-1)*`util'^(`tet'-1))/`u')
			local d1copd1v		(`cop'*(`A'^(`tettil'-1)*`vtil'^(`tet'-1))/`v')
			local d1copd1tet	(`tettil'*`A'^(`tettil')*`cop'*(`tettil'*ln(`A') - (`d1Ad1tet'/`A')))	
		}
		if "`copfc'" == "joe"{
			local d1copd1u 		(`A'^(`tettil'-1)*`util'^(`tet'-1)*(1-`vtil'^(`tet')))
			local d1copd1v		(`A'^(`tettil'-1)*`vtil'^(`tet'-1)*(1-`util'^(`tet')))
			local d1copd1tet	(`tettil'^2*`A'^(`tettil')*ln(`A') - `tettil'*`A'^(`tettil'-1)*`d1Ad1tet')
		}
		if "`copfc'" == "amh"{
			local d1copd1u 		((`v'/`A') - `u'*`v'*`tet'*(`vtil'/(`A'^2)))
			local d1copd1v		((`u'/`A') - `u'*`v'*`tet'*(`util'/(`A'^2)))
			local d1copd1tet	(`u'*`v'*`util'*(`vtil'/(`A'^2)))	
		}	
		
		
		*! second derivatives of copulas
		if "`copfc'" == "product"{
			local d2copd2u		(0)	
			local d2copd1ud1v	(1)
			local d2copd1ud1tet	(0)
			
			local d2copd2v		(0)
			local d2copd1vd1tet	(0)
			
			local d2copd2tet	(0)
		}
		if "`copfc'" == "gaussian"{
			
		}
		if "`copfc'" == "fgm"{
			local d2copd2u		(-2*`tet'*`vtil'*`v')	
			local d2copd1ud1v	(1 + `tet'*(1-2*`v')*(1-2*`u'))
			local d2copd1ud1tet	(`vtil'*`v'*(1-2*`u'))
			
			local d2copd2v		(-2*`tet'*`util'*`u')	
			local d2copd1vd1tet	(`util'*`u'*(1-2*`v'))
			
			local d2copd2tet	(0)		
		}
		if "`copfc'" == "plackett"{
			local d2copd2u		(.5*((1-`tet')/sqrt(`A') + ///
								(`tettil'*(`r'-2*`v'*`tet')^2)/(`A'*sqrt(`A'))))	
			local d2copd1ud1v	(.5*((1+`tet')/sqrt(`A') + ///
								(`tettil'*(`r'-2*`v'*`tet')*(`r'-2*`u'*`tet'))/(`A'*sqrt(`A'))))
			local d2copd1ud1tet	(.5*( (`v'-`u')/sqrt(`A') - ///
								(2*`v'*`tet' - `r')/(2*`A'*sqrt(`A')) * `d1Ad1tet'))
			
			local d2copd2v		(.5*((1-`tet')/sqrt(`A') + ///
								(`tettil'*(`r'-2*`u'*`tet')^2)/(`A'*sqrt(`A'))))	
			local d2copd1vd1tet	(.5*( (`u'-`v')/sqrt(`A') - ///
								(2*`u'*`tet' - `r')/(2*`A'*sqrt(`A')) * `d1Ad1tet'))
			
			local d2copd2tet	(`d1copd1tet'*(1/(`tettil'*sqrt(`A')) - 1/`tettil' - (`d1Ad1tet')/(2*`A')))
		}
		if "`copfc'" == "clayton"{
			local d2copd2u		(`d1copd1u'*((1+`tet')/(`A'*`u'^(`tet'+1)) - (1+`tet')/(`u'))) 
			local d2copd1ud1v	((1+`tet')/(`A'^(`tettil'+2)*(`u'*`v')^(`tet'+1)))
			local d2copd1ud1tet	(`d1copd1u'*(((-`tettil'-1)/`A')*`d1Ad1tet' ///
								+ ln(`A')*(`tettil'^2) - ln(`u')))
			
			local d2copd2v		(`d1copd1v'*((1+`tet')/(`A'*`v'^(`tet'+1)) - (1+`tet')/(`v'))) 
			local d2copd1vd1tet	(`d1copd1v'*(((-`tettil'-1)/`A')*`d1Ad1tet' ///
								+ ln(`A')*(`tettil'^2) - ln(`v')))
			
			local d2copd2tet	((`d1copd1tet'^2)/`cop' - `d1copd1tet'*`tettil' + `cop'*`tettil'* ///
								((`tettil'/`A')*`d1Ad1tet' - ln(`A')*`tettil'^2 - ///
								(`d2Ad2tet'/`A') + ((`d1Ad1tet'^2)/(`A'^2))))
		}
		if "`copfc'" == "frank"{
			local d2copd2u		(`tet'*(`d1copd1u'^2) - `tet'*`d1copd1u')	
			local d2copd1ud1v	((-`tet'*`util'*`vtil')/(`A'*`A'*(`eta'-1)))
			local d2copd1ud1tet	(`d1copd1u'*(`eta'/(`eta'-1) - (`d1Ad1tet'/`A') - `u' - (`vtil'*`v')/(`vtil'-1)))
			
			local d2copd2v		(`tet'*(`d1copd1v'^2) - `tet'*`d1copd1v')	
			local d2copd1vd1tet	(`d1copd1v'*(`eta'/(`eta'-1) - (`d1Ad1tet'/`A') - `v' - (`util'*`u')/(`util'-1)))
			
			local d2copd2tet	((`tettil'/`A')*`d1Ad1tet'*(2*`tettil' + (`d1Ad1tet'/`A')) + ///
								`tettil'*((`d2Ad2tet'/`A') - 2*(`tettil'^2)*ln(`A')))
		}
		if "`copfc'" == "gumbel"{
			local d2copd2u		(((1-`tet')/(`util'*`u'))*`d1copd1u'*(1 - (`util'^(`tet')/`A')) ///
								+ ((`d1copd1u'^2)/`cop') - (`d1copd1u'/`u'))	
			local d2copd1ud1v	((`util'^(`tet'-1)/(`A'*`u'))*`d1copd1v'*(`A'^(`tettil') - 1 + `tet'))
			local d2copd1ud1tet	((`d1copd1u'*`d1copd1tet'/`cop')*(1 - `A'^(-`tettil')) ///
								- `d1copd1u'*((`d1Ad1tet'/`A') - ln(`util')))
		
			local d2copd2v		(((1-`tet')/(`vtil'*`v'))*`d1copd1v'*(1 - (`vtil'^(`tet')/`A')) ///
								+ (`d1copd1v'^2/`cop') - (`d1copd1v'/`v'))	
			local d2copd1vd1tet	((`d1copd1v'*`d1copd1tet'/`cop')*(1 - `A'^(-`tettil')) ///
								- `d1copd1v'*((`d1Ad1tet'/`A') - ln(`vtil')))
			
			local d2copd2tet	((`d1copd1tet'^2/`cop')*(1 - `A'^(-`tettil')) ///
								- 2*`tettil'*`d1copd1tet' + (`tettil'*`A'^(`tettil'-1)*`cop')*( ///
								(`d1Ad1tet'^2/`A') - `d2Ad2tet'))
		}
		if "`copfc'" == "joe"{
			local d2copd2u		(((1-`tet')/ `util')*`d1copd1u' - ((1-`tet')/`A'^(`tettil'))*`d1copd1u'^2)	
			local d2copd1ud1v	(`A'^(`tettil'-2)*(`util'*`vtil')^(`tet'-1)*(`A'-1+`tet'))
			local d2copd1ud1tet	(`d1copd1u'*(ln(`util') - (`d1Ad1tet'/`A') ///
								- ((`vtil'^(`tet')*ln(`vtil'))/(1-`vtil'^(`tet')))) ///
								- (`d1copd1u'*`d1copd1tet'/`A'^(`tettil')))
			
			local d2copd2v		(((1-`tet')/ `vtil')*`d1copd1v' - ((1-`tet')/`A'^(`tettil'))*`d1copd1v'^2)	
			local d2copd1vd1tet	(`d1copd1v'*(ln(`vtil') - (`d1Ad1tet'/`A') ///
								- ((`util'^(`tet')*ln(`util'))/(1-`util'^(`tet')))) ///
								- (`d1copd1v'*`d1copd1tet'/`A'^(`tettil')))
			
			local d2copd2tet	(`tettil'*`A'^(`tettil'-1)*((`d1Ad1tet'^2/`A') - `d2Ad2tet') ///
								- 2*`tettil'*`d1copd1tet' - (`d1copd1tet'^2/`A'^(`tettil')))
		}
		if "`copfc'" == "amh"{
			local d2copd2u		(-2*(`tet'*`vtil'/`A')*`d1copd1u')	
			local d2copd1ud1v	((1/`A') + (`u'*`v'*`tet'/`A'^2) - (`d1copd1u'*`tet'*`util'/`A') ///
								- (`d1copd1v'*`tet'*`vtil'/`A'))
			local d2copd1ud1tet	(`d1copd1tet'*((1/`u') - (1/`util') - (2*`tet'*`vtil'/`A')))
			
			local d2copd2v		(-2*(`tet'*`util'/`A')*`d1copd1v')	
			local d2copd1vd1tet	(`d1copd1tet'*((1/`v') - (1/`vtil') - (2*`tet'*`util'/`A')))
			
			local d2copd2tet	(`d1copd1tet'*(2*`util'*`vtil'/`A'))
		}	
	}
	
	
	
	if `"`d1'`d2'"' == ""{
			
		*!	joint probabilities
		if inlist("`type'","","p11","p10","p01","p00"){
			if "`type'" == "" {
				local type "p11"
				di in gr "(option p11 assumed; Pr(`dep1'=1,`dep2'=1))"
			}	
		
			local val1 = substr("`type'",2,1)
			local val2 = substr("`type'",3,1)
			local pred "Pr(`dep1'=`val1',`dep2'=`val2')"
		
			gen `vtyp' `varn' = `s'*`t' - `t'*`q1'*`u' - `s'*`q2'*`v' + `q1'*`q2'*`cop'		`if' `in'
			label var `varn' "`pred'"
			exit
		}
		
						
		*!	marginal probability of xb1
		if "`type'" == "pmarg1"{
			tempvar	xb
			
			_predict double `xb' `if' `in', eq(#1) `offset'
			
			gen `vtyp'	`varn' = normal(`xb')	`if' `in'
			label var	`varn' "Pr(`dep1'=1)"
			exit
		}	
		
		
		*!	marginal probability of xb2
		if "`type'" == "pmarg2"{
			tempvar	zg
			
			_predict double `zg' `if' `in', eq(#2) `offset'
			
			gen `vtyp'	`varn' = normal(`zg')	`if' `in'
			label var	`varn' "Pr(`dep2'=1)"
			exit
		}		
	
	
		*!	conditional probabilities
		if inlist("`type'","pcond1","pcond10"){
			if "`type'" == "pcond1"{
				local pred	"Pr(`dep1'=1|`dep2'=1)"
			}
			else{
				local pred	"Pr(`dep1'=1|`dep2'=0)"
			}
			
			gen `vtyp'	`varn' = 	(`s'*`t' - `t'*`q1'*`u' - `s'*`q2'*`v' ///
									+ `q1'*`q2'*`cop') / normal(`q2'*`zg')	`if' `in'
			label var	`varn' "`pred'"
			exit
		}	
	
	
		*!	conditional probability for atec estimation
		if "`type'" == "pcond2"{

			gen `vtyp'	`varn' = (`s'*`t' - `t'*`q1'*`u' - `s'*`q2'*`v' ///
									+ `q1'*`q2'*`cop') / normal(`q1'*`xb')	`if' `in'
			label var	`varn' "Pr(`dep2'=1|`dep1'=1)"			
			exit
		}	
	
	
			
		*!	special joint probability for atet estimation
		if "`type'" == "pmargcond1"{
				
			gen `vtyp'	`varn' = (`s'*`t' - `t'*`q1'*`u' - `s'*`q2'*`v' ///
									+ `q1'*`q2'*`cop')	`if' `in'
			label var	`varn' "Pr(`dep1'=1|`dep2'=1): Conditional Marginal Probability"
			exit
		}	
		
		
	}
	
	
	
	if `"`d1'`d2'"' != ""{
	
		*!	auxiliary derivatives
		local d1Prd1u		(-`t'*`q1' + `q1'*`q2'*`d1copd1u')
		local d1Prd1v		(-`s'*`q2' + `q1'*`q2'*`d1copd1v')
		local d1Prd1tet		(`q1'*`q2'*`d1copd1tet')
		
		local d2Prd2u		(`q1'*`q2'*`d2copd2u')
		local d2Prd1ud1v	(`q1'*`q2'*`d2copd1ud1v')
		local d2Prd1ud1tet	(`q1'*`q2'*`d2copd1ud1tet')
		
		local d2Prd2v		(`q1'*`q2'*`d2copd2v')
		local d2Prd1vd1tet	(`q1'*`q2'*`d2copd1vd1tet')

		local d2Prd2tet		(`q1'*`q2'*`d2copd2tet')
	
		local w1			(`q1'*`xb')
		local w2			(`q2'*`zg')
		local Pr1			(normal(`w1'))
		local Pr2			(normal(`w2'))
		local BiPr			(`s'*`t' - `t'*`q1'*`u' - `s'*`q2'*`v' + `q1'*`q2'*`cop')
		
		local d1Pr1d1xb		(`q1'*normalden(`w1'))
		local d1Pr2d1zg		(`q2'*normalden(`w2'))
		local d2Pr1d2xb		(-`w1'*normalden(`w1'))
		local d2Pr2d2zg		(-`w2'*normalden(`w2'))
		
		
				
		*!	joint probability
		if inlist("`type'","","p11","p10","p01","p00","pmargcond1"){
			
			if "`type'" == "" {
				local type "p11"
				local val1 = substr("`type'",2,1)
				local val2 = substr("`type'",3,1)
				local pred "Pr(`dep1'=`val1',`dep2'=`val2')"
			}	
			else if "`type'" == "pmargcond1"{
				local pred "Pr(`dep1'=1|`dep2'=1): Conditional Marginal Probability"
			}
			else{
				local val1 = substr("`type'",2,1)
				local val2 = substr("`type'",3,1)
				local pred "Pr(`dep1'=`val1',`dep2'=`val2')"
			}
			
			
			*!	first derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & `"`d2'"' == ""{
				if `"`d1'"' == "#1"{
					gen `vtyp' 	`varn' = `d1Prd1u'*`d1ud1xb'		`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#2"{
					gen `vtyp' 	`varn' = `d1Prd1v'*`d1vd1zg'		`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#3"{
					gen `vtyp'	`varn' = `d1Prd1tet'*`d1tetd1del'	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
			}
			
			
			*!	second derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & inlist(`"`d2'"',"#1","#2","#3"){
				if `"`d1'`d2'"' == "#1#1"{
					gen `vtyp' 	`varn' =  `d2Prd2u'*(`d1ud1xb'^2) + `d1Prd1u'*`d2ud2xb'		`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#2", "#2#1"){
					gen `vtyp'	`varn' =  `d2Prd1ud1v'*`d1vd1zg'*`d1ud1xb'		`if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#3", "#3#1"){
					gen `vtyp' 	`varn' = `d2Prd1ud1tet'*`d1tetd1del'*`d1ud1xb'	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#2#2"{
					gen `vtyp'	`varn' = `d2Prd2v'*(`d1vd1zg'^2) + `d1Prd1v'*`d2vd2zg'		`if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#2#3", "#3#2"){
					gen `vtyp' 	`varn' = `d2Prd1vd1tet'*`d1tetd1del'*`d1vd1zg'	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#3#3"{
					gen `vtyp'	`varn' = `d2Prd2tet'*(`d1tetd1del'^2) + `d1Prd1tet'*`d2tetd2del' 	`if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
			}		
			exit		
		}
		
		
		
		*!	marginal probability of xb1
		if "`type'" == "pmarg1"{
			tempvar	xb
			_predict double `xb' `if' `in', eq(#1) `offset'
			local pred "Pr(`dep1'=1)"		

			
			*!	first derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & `"`d2'"' == ""{
				if `"`d1'"' == "#1"{
					gen `vtyp' 	`varn' = normalden(`xb')	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#2"{
					gen `vtyp' 	`varn' = 0	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#3"{
					gen `vtyp'	`varn' = 0	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
			}
			
			
			*!	second derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & inlist(`"`d2'"',"#1","#2","#3"){
				if `"`d1'`d2'"' == "#1#1"{
					gen `vtyp' 	`varn' =  -`xb'*normalden(`xb')		`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#2", "#2#1"){
					gen `vtyp'	`varn' = 0	`if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#3", "#3#1"){
					gen `vtyp' 	`varn' = 0	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#2#2"{
					gen `vtyp'	`varn' = 0	`if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#2#3", "#3#2"){
					gen `vtyp' 	`varn' = 0	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#3#3"{
					gen `vtyp'	`varn' = 0 `if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
			}
			exit
		}
		

		*!	marginal probability of xb2
		if "`type'" == "pmarg2"{
			tempvar	zg
			_predict double `zg' `if' `in', eq(#2) `offset'
			local pred "Pr(`dep2'=1)"
		
		
			*!	first derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & `"`d2'"' == ""{
				if `"`d1'"' == "#1"{
					gen `vtyp' 	`varn' = 0	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#2"{
					gen `vtyp' 	`varn' = normalden(`zg')	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#3"{
					gen `vtyp'	`varn' = 0		`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
			}
			
			
			*!	second derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & inlist(`"`d2'"',"#1","#2","#3"){
				if `"`d1'`d2'"' == "#1#1"{
					gen `vtyp' 	`varn' = 0	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#2", "#2#1"){
					gen `vtyp'	`varn' = 0	`if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#3", "#3#1"){
					gen `vtyp' 	`varn' = 0	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#2#2"{
					gen `vtyp'	`varn' = -`zg'*normalden(`zg')	`if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#2#3", "#3#2"){
					gen `vtyp' 	`varn' = 0	`if' `in'
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#3#3"{
					gen `vtyp'	`varn' = 0 `if' `in'
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
			}
			exit
		}
		

		*!	conditional probabilities
		if inlist("`type'","pcond1","`pcond10'"){
			if "`type'" == "pcond1"{
				local pred "Pr(`dep1'=1|`dep2'=1)"
			}
			else{
				local pred "Pr(`dep1'=1|`dep2'=0)"
			}
		
		
			*!	first derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & `"`d2'"' == ""{
				if `"`d1'"' == "#1"{
					gen `vtyp' 	`varn' = (1/`Pr2')*(`d1Prd1u'*`d1ud1xb')	`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#2"{
					gen `vtyp' 	`varn' = (1/(`Pr2'^2)) * (`Pr2'*(`d1Prd1v'*`d1vd1zg') - ///
										`BiPr'*`d1Pr2d1zg') 	`if' `in'
										
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#3"{
					gen `vtyp'	`varn' = (1/`Pr2')*(`d1Prd1tet'*`d1tetd1del')	`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
			}
			
			
			*!	second derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & inlist(`"`d2'"',"#1","#2","#3"){
				if `"`d1'`d2'"' == "#1#1"{
					gen `vtyp' 	`varn' = (1/`Pr2')*(`d2Prd2u'*(`d1ud1xb'^2) + ///
										`d1Prd1u'*`d2ud2xb') 	`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#2", "#2#1"){
					gen `vtyp'	`varn' = (1/(`Pr2'^2))*(`Pr2'*(`d2Prd1ud1v'*`d1vd1zg'*`d1ud1xb') - ///
										(`d1Prd1u'*`d1ud1xb')*`d1Pr2d1zg') 	`if' `in'
					
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#3", "#3#1"){
					gen `vtyp' 	`varn' = (1/`Pr2')*(`d2Prd1ud1tet'*`d1tetd1del'*`d1ud1xb')	`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#2#2"{
					gen `vtyp'	`varn' = (1/(`Pr2'^2))*( ///
								`Pr2'*(`d2Prd2v'*(`d1vd1zg'^2) + `d1Prd1v'*`d2vd2zg') - ///
								2*(`d1Prd1v'*`d1vd1zg')*`d1Pr2d1zg' + ///
								2*(`BiPr'/`Pr2')*(`d1Pr2d1zg'^2) - `BiPr'*`d2Pr2d2zg') ///
								`if' `in'						 
					
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#2#3", "#3#2"){
					gen `vtyp' 	`varn' = (1/(`Pr2'^2))*(`Pr2'*(`d2Prd1vd1tet'*`d1tetd1del'*`d1vd1zg') - ///
										(`d1Prd1tet'*`d1tetd1del')*`d1Pr2d1zg') 	`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#3#3"{
					gen `vtyp'	`varn' = (1/`Pr2')*(`d2Prd2tet'*(`d1tetd1del'^2) + ///
								`d1Prd1tet'*`d2tetd2del')	`if' `in'
					
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
			}
			exit
		}
		
		
		*!	conditional probability fpr atec estimation
		if "`type'" == "pcond2"{
			local pred "Pr(`dep2'=1|`dep1'=1)"
		
		
			*!	first derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & `"`d2'"' == ""{
				if `"`d1'"' == "#1"{
					gen `vtyp' 	`varn' = (1/(`Pr1'^2))*( `Pr1'*(`d1Prd1u'*`d1ud1xb') - ///
								`BiPr'*(`d1Pr1d1xb')) 	`if' `in'
								
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#2"{
					gen `vtyp' 	`varn' = (1/`Pr1')*(`d1Prd1v'*`d1vd1zg')		`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
				if `"`d1'"' == "#3"{
					gen `vtyp'	`varn' = (1/`Pr1')*(`d1Prd1tet'*`d1tetd1del')	`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1')"
				}
			}
			
			
			*!	second derivatives
			if inlist(`"`d1'"',"#1","#2","#3") & inlist(`"`d2'"',"#1","#2","#3"){
				if `"`d1'`d2'"' == "#1#1"{
					gen `vtyp' 	`varn' = (1/(`Pr1'^2))*( ///
								`Pr1'*(`d2Prd2u'*(`d1ud1xb'^2) + `d1Prd1u'*`d2ud2xb') - ///
								2*(`d1Prd1u'*`d1ud1xb')*`d1Pr1dxb' + ///
								2*(`BiPr'/`Pr1')*(`d1Pr1d1xb'^2) - `BiPr'*`d2Pr1d2xb') ///
								`if' `in'
								
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#2", "#2#1"){
					gen `vtyp'	`varn' = (1/(`Pr1'^2))*(`Pr1'*(`d2Prd1ud1v'*`d1vd1zg'*`d1ud1xb') -  ///
								`d1Pr1d1xb'*(`d1Prd1v'*`d1vd1zg')		`if' `in'
								
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#1#3", "#3#1"){
					gen `vtyp' 	`varn' = (1/(`Pr1'^2))*(`Pr1'*(`d2Prd1ud1tet'*`d1tetd1del'*`d1ud1xb') -  ///
								`d1Pr1d1xb'*(`d1Prd1tet'*`d1tetd1del')	`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#2#2"{
					gen `vtyp'	`varn' = (1/`Pr1')*(`d2Prd2v'*(`d1vd1zg'^2) + ///
								`d1Prd1v'*`d2vd2zg')	`if' `in'
					
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if inlist(`"`d1'`d2'"', "#2#3", "#3#2"){
					gen `vtyp' 	`varn' = (1/`Pr1')*(`d2Prd1vd1tet'*`d1tetd1del'*`d1vd1zg')	`if' `in'
					
					label var 	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
				if `"`d1'`d2'"' == "#3#3"{
					gen `vtyp'	`varn' = (1/`Pr1')*(`d2Prd2tet'*(`d1tetd1del'^2) + ///
								`d1Prd1tet'*`d2tetd2del')	`if' `in'
					
					label var	`varn' "d `pred' / d xb(`d1') d xb(`d2')"
				}
			}
			exit
		}
	}
	
	error 198			//	if user specified more than one option	
end


program define rmTS, rclass
	
	version 11
	
	local tsnm = cond( match("`0'", "*.*"),  		/*
			*/ bsubstr("`0'", 			/*
			*/	  (index("`0'",".")+1),.),     	/*
			*/ "`0'")

	return local rmTS `tsnm'
end
