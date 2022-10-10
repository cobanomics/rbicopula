*! version 1.1.0 , 10oct2022
*! Author: Mustafa Coban, Institute for Employment Research (Germany)
*! Website: mustafacoban.de
*! Support: mustafa.coban@iab.de


/****************************************************************/
/*    			 rbicopula lf1 evaluator						*/
/****************************************************************/


program define rbicopula_lf1

	version 11
	args todo b lnfi g1 g2 g3 H
	
	tempvar 	xb zg
	tempname 	delta
	
	local copfunc = "${COPULAFUNC}"
	
	
	mleval `xb'  = `b', eq(1)
	mleval `zg'  = `b', eq(2)
	mleval `delta' = `b', eq(3) scalar
	

	
	*!	common notational shortcuts	
	tempvar  s t q1 q2 u v d1ud1xb d1vd1zg
	
	qui{
		gen byte `s' = $ML_y1
		gen byte `t' = $ML_y2
		gen byte `q1' = 2*$ML_y1 - 1
		gen byte `q2' = 2*$ML_y2 - 1		
		
		gen double `u' = normal(-`xb')
		gen double `v' = normal(-`zg')
		
		gen double `d1ud1xb' = -normalden(-`xb')
		gen double `d1vd1zg' = -normalden(-`zg')
	}

	
	*!	dependence parameter related
	tempname	theta d1td1d d2td2d
	
	if inlist("`copfunc'","product"){
		scalar `theta' 	= 0
		scalar `d1td1d' = 0
	}
	else if inlist("`copfunc'","gaussian","fgm","amh"){
		scalar `theta' = (exp(2*`delta')-1) / (1+exp(2*`delta'))
		scalar `d1td1d' = (4 * exp(2*`delta')) / ( (1+exp(2*`delta'))*(1+exp(2*`delta')) )
	}
	else if inlist("`copfunc'","plackett","clayton"){
		scalar `theta' = exp(`delta')
		scalar `d1td1d' = exp(`delta')
	}	
	else if inlist("`copfunc'","frank"){
		scalar `theta' = `delta'
		scalar `d1td1d' = 1
	}
	else if inlist("`copfunc'","gumbel","joe"){
		scalar `theta' = (1 + exp(`delta'))
		scalar `d1td1d' = exp(`delta')
	}
	
	
	
	*!	notational shortcuts for copula
	qui{
		if "`copfunc'" == "gaussian"{
			tempname	eta
			tempvar		w1 w2
			
			scalar `eta' = 1/(sqrt(1-`theta'*`theta'))
			gen double `w1' = invnormal(`u')
			gen double `w2' = invnormal(`v')
		}
		if "`copfunc'" == "fgm"{
			tempvar		util vtil
			
			gen double `util' = (1 - `u')
			gen double `vtil' = (1 - `v')
		}
		if "`copfunc'" == "plackett"{
			tempvar		r A
			tempname	thetatil
			
			scalar `thetatil' 	= (`theta' - 1)
			gen double `r' 		= 1 + `thetatil'*(`u'+`v')
			gen double `A' 		= `r'*`r' - 4*`u'*`v'*`theta'*`thetatil'
		}
		if "`copfunc'" == "clayton"{
			tempvar		A
			
			gen double `A' = (`u'^(-`theta') + `v'^(-`theta') - 1)
		}
		if "`copfunc'" == "frank"{
			tempvar		util vtil A
			tempname	eta
			
			scalar `eta' 		= (exp(-`theta') - 1)
			gen double `util' 	= (exp(-`theta'*`u') - 1)
			gen double `vtil' 	= (exp(-`theta'*`v') - 1)
			gen double `A' 		= 1 + ((`util'*`vtil')/(`eta'))
		}
		if "`copfunc'" == "gumbel"{
			tempvar		util vtil A
			
			gen double `util' 	= -ln(`u')
			gen double `vtil' 	= -ln(`v')
			gen double `A' 		= (`util'^`theta' + `vtil'^`theta')
		}
		if "`copfunc'" == "joe"{
			tempvar		util vtil A
			tempname	eta
			
			scalar `eta' 		= 1/`theta'
			gen double `util' 	= (1 - `u')
			gen double `vtil' 	= (1 - `v')
			gen double `A' 		= (`util'^`theta' + `vtil'^`theta' - (`util'^`theta')*(`vtil'^`theta'))
		}
		if "`copfunc'" == "amh"{
			tempvar		util vtil A
			
			gen double `util' 	= (1 - `u')
			gen double `vtil' 	= (1 - `v')
			gen double `A' 		= (1 - `theta'*`util'*`vtil')
		}
	}
	
	
	
	*!	copula functions
	tempvar		cop
	
	qui{
		if "`copfunc'" == "product"{
			gen double `cop' = `u'*`v'
		}
		if "`copfunc'" == "gaussian"{
			gen double `cop' = binormal(`w1', `w2', `theta')
		}
		if "`copfunc'" == "fgm"{
			gen double `cop' = `u'*`v'*(1 + `theta'*`util'*`vtil')
		}
		if "`copfunc'" == "plackett"{
			gen double `cop' = (`r' - sqrt(`A')) / (2*`thetatil')
		}
		if "`copfunc'" == "clayton"{
			gen double `cop' = `A'^(-1/`theta')
		}
		if "`copfunc'" == "frank"{
			gen double `cop' = (-1/`theta') * ln(`A')
		}
		if "`copfunc'" == "gumbel"{
			gen double `cop' = exp(-(`A'^(1/`theta')))
		}
		if "`copfunc'" == "joe"{
			gen double `cop' = 1 - (`A'^(`eta'))
		}
		if "`copfunc'" == "amh"{
			gen double `cop' = (`u'*`v')/`A'
		}
	}
	
	
	tempvar		prob ell
	
	qui{
		gen double `prob'	= `s'*`t' - `t'*`q1'*`u' - `s'*`q2'*`v' + `q1'*`q2'*`cop'
		gen double `ell'	= `prob'
				
		replace `lnfi'		= ln(`ell')
	}
	
	if (`todo' == 0) exit			//	end of lf0-evaluator
	
	
	
	
	*!	first derivatives of copula functions
	tempvar		d1copd1u d1copd1v d1copd1t
	
	qui{
		if "`copfunc'" == "product"{
			gen double `d1copd1u' = `v'
			gen double `d1copd1v' = `u'
			gen double `d1copd1t' = 0
		}
		if "`copfunc'" == "gaussian"{
			gen double `d1copd1u' = normal(`eta'*(`w2' - `theta'*`w1') )
			gen double `d1copd1v' = normal(`eta'*(`w1' - `theta'*`w2') )
			gen double `d1copd1t' = (1/(2*_pi))*`eta'* ///
						exp(-.5*(`eta')^2*((`w1'^2) - 2*`theta'*`w1'*`w2' + (`w2'^2)))
		}
		if "`copfunc'" == "fgm"{
			gen double `d1copd1u' = `v'*(1+`theta'*`util'*`vtil') - `theta'*`u'*`v'*`vtil'
			gen double `d1copd1v' = `u'*(1+`theta'*`util'*`vtil') - `theta'*`u'*`v'*`util'
			gen double `d1copd1t' = `u'*`v'*`util'*`vtil'
		}
		if "`copfunc'" == "plackett"{
			gen double `d1copd1u' = .5 - ( (`r' - 2*`v'*`theta')/(2*sqrt(`A')) )
			gen double `d1copd1v' = .5 - ( (`r' - 2*`u'*`theta')/(2*sqrt(`A')) )
			gen double `d1copd1t' = (1/(`thetatil'*sqrt(`A'))) * (`cop' - `u'*`v')
		}
		if "`copfunc'" == "clayton"{
			gen double `d1copd1u' = `A'^((-1/`theta') - 1) * `u'^(-`theta'-1)
			gen double `d1copd1v' = `A'^((-1/`theta') - 1) * `v'^(-`theta'-1)
			gen double `d1copd1t' = `cop'*( (1/(`A'*`theta')) * ( ///
						`u'^(-`theta')*ln(`u') + `v'^(-`theta')*ln(`v')) ///
						+ (ln(`A')/(`theta'*`theta')) )
		}
		if "`copfunc'" == "frank"{
			gen double `d1copd1u' = (`vtil'*exp(-`theta'*`u')) / (`A'*`eta') 
			gen double `d1copd1v' = (`util'*exp(-`theta'*`v')) / (`A'*`eta') 
			gen double `d1copd1t' = (ln(`A')/(`theta'*`theta'))  - ///
						(1/(`theta'*`A'*`eta'*`eta')) * ( ///
						(-`u'*`vtil'*exp(-`theta'*`u') - `v'*`util'*exp(-`theta'*`v'))*`eta' ///
						+ `util'*`vtil'*exp(-`theta') )
		}
		if "`copfunc'" == "gumbel"{
			gen double `d1copd1u' = `cop'*((`util'^(`theta'-1)) / `u') * `A'^((1/`theta') - 1)
			gen double `d1copd1v' = `cop'*((`vtil'^(`theta'-1)) / `v') * `A'^((1/`theta') - 1)
			gen double `d1copd1t' = `cop'* (-(`A'^(1/`theta'))/`theta') * ( ///
						((`util'^(`theta')*ln(`util') + `vtil'^(`theta')*ln(`vtil')) /`A') ///
						- (ln(`A')/`theta') )
		}
		if "`copfunc'" == "joe"{
			gen double `d1copd1u' = `A'^(`eta'-1) * `util'^(`theta'-1) * (1 - `vtil'^(`theta'))
			gen double `d1copd1v' = `A'^(`eta'-1) * `vtil'^(`theta'-1) * (1 - `util'^(`theta'))
			gen double `d1copd1t' = -`eta'*`A'^(`eta') * ( ///
						(1/`A') * (`util'^(`theta')*ln(`util')*(1-`vtil'^(`theta')) + ///
							`vtil'^(`theta')*ln(`vtil')*(1-`util'^(`theta'))) ///
						- ln(`A')*`eta')
		}
		if "`copfunc'" == "amh"{
			gen double `d1copd1u' = (`v'/`A') * (1 - ((`u'*`theta'*`vtil')/`A'))
			gen double `d1copd1v' = (`u'/`A') * (1 - ((`v'*`theta'*`util')/`A'))
			gen double `d1copd1t' = (`u'*`v'*`util'*`vtil')/(`A'*`A')
		}
	}	
	
	
	
	*!	scores
	tempvar		scr1 scr2 scr3
	
	qui{
		gen double `scr1' = (1/`ell') * (-`t'*`q1' + `q1'*`q2'*`d1copd1u') * `d1ud1xb'
		gen double `scr2' = (1/`ell') * (-`s'*`q2' + `q1'*`q2'*`d1copd1v') * `d1vd1zg' 
		gen double `scr3' = (1/`ell') * (`q1'*`q2'*`d1copd1t') * `d1td1d' 	
		
		replace `g1' = `scr1'
		replace `g2' = `scr2'
		replace `g3' = `scr3'
	}
	
	if (`todo' == 1) exit			//	end of lf1-evaluator
	
end
