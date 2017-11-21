## -------------------------------------------------------------------------- ##
##
## THIS IS THE CONTROL FILE
## 
## -------------------------------------------------------------------------- ##

## -------------------------------------------------------------------------- ##
## This file controls estimation phases, model years, and blocks for changes 
## to fecundity, maturity-at-age, natural mortality, and gear selectivity.
## -------------------------------------------------------------------------- ##

#Number of age classes (nages)*
6
#First year of data (dat_styr)*
1974
#Last year of data (dat_endyr)*
2015
#Start year for modeling (mod_styr)*
1988
#Last year of model (mod_endyr)*
2015

#Current Harvest Threshold (tons) 
#
5000
#Threhold (tons)
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000
5000



#Weight-at-age used for biomass forecast
#
59.59707792	70.08571429	88.73488372	101.732	107.0680851	122.5914286


#
#Year blocks in Fecundity*
#
1
#					
#Change in Fecundity years*
#
1988
2015
#					
#Fecundity slope*
#
210.5
#					
#Fecundity intercept*
#
1092.3


#
#Define Parameter structures
#Year blocks in survival matrix (S_Bk)
3
#Years for survival breaks (S_Bk_Yrs)*
1988
1999
2008
2015
#Year blocks in maturity matrix (mat_Bk)*
3
#Years for maturity breaks (mat_Bk_Yrs)*
1988
1999
2003
2015
#Year blocks in gear selectivity matrix (gs_Bk)*
1
#Years for selectivity breaks (gs_Bk_Yrs)*
1988
2015
#
#Parameter phases
#Initial population (ph_Int)*
1
#Survival/Mortality (ph_S)
1
#Maturity inflection age (ph_mat_a)
1
#Gear selectivity inflection age (ph_gs_a)
1
#Gear selectivity slope (ph_gs_b)
1
#Maturity slope (ph_mat_b)
2
#Recruitment (ph_Rec)
3
#Ricker spawner-recruit (ph_Ric)
3

#Mile-days coefficient
-3
#
#Objective function weighting
#Catch age composition*(lC)
1
#Spawning age composition*(lS)
1
#Ricker spawner-recruit*(lR)#recheck every few years to make sure no influence on model output (or weight less)
0.001

#Egg deposition*(lE)(fill in new weights)
0
0
0
0
0
0
0
0
0
0
0
0
0
0
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
#Mile-days of Milt*
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0

#EOF
42


