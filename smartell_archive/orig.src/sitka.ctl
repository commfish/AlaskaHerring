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
1971
#Last year of data (dat_endyr)*
2015
#Start year for modeling (mod_styr)*
1980
#Last year of model (mod_endyr)*
2015

#Harvest Threshold (tons) Thresh //changed to tons
#
25000
#Threshold calculation denominator (tons)//changed to tons
#
20000
#Threshold vector (tons)
7500
7500
7500
7500
7500
7500
7500
7500
7500
7500
7500
7500
7500
7500
7500
7500
7500
20000
20000
20000
20000
20000
20000
20000
20000
20000
20000
20000
20000
20000
25000
25000
25000
25000
25000
25000
#
#Forecast weight-at-age (fw_a_a)
#
74.85279503	88.94545455	113.4243697	141.112	158.4052632	182.4891089

#Year blocks in Fecundity*
#
4
#					
#Change in Fecundity years*
#
1980
1984
1998
2002
2015
#					
#Fecundity slope*
#
195.314028
222.55
226.1260221
198.78295
#					
#Fecundity intercept*
#
1840.107419
2503.9
3383.01566
3031.97919
#

#
#Define Parameter structures
#Year blocks in survival matrix (S_Bk)*
2
#Years for survival breaks (S_Bk_Yrs)*
1980
1999
2015
#Year blocks in maturity matrix (mat_Bk)*
1
#Years for maturity breaks (mat_Bk_Yrs)*
1980
2015
#Year blocks in gear selectivity matrix (gs_Bk)*
1
#Years for selectivity breaks (gs_Bk_Yrs)*
1980
2015
#
#  

#Parameter phases
#Initial population (ph_Int)*
1
#Survival/Mortality (ph_S)
2
#Maturity inflection age (ph_mat_a)
1
#Gear selectivity inflection age (ph_gs_a)
1
#Gear selectivity slope (ph_gs_b)
2
#Maturity slope (ph_mat_b)
3
#Recruitment (ph_Rec)
4
#Ricker spawner-recruit (ph_Ric)
3
#Mile-days coefficient (ph_md)
1
#
#Objective function weighting
#Catch*
#1
#Catch age composition*(lC)
1
#Spawning age composition*(lS)
1
#Ricker spawner-recruit*(lR)
0.001

#Egg deposition* (lE)
0.75
0.75
0.75
0.75
0.75
0.75
0.75
607.1983401
8.004756537
2.838738274
3.193464659
3.266557274
4.387542317
2.632507984
3.130142972
4.595565238
1.416959688
0.631767624
4.118034845
6.239086089
3.17750044
1.280179608
1.577145861
13.44689397
1.802620908
1.884634303
1.797786974
2.438931415
0.803652368
0.549103353
0.925326911
1.727719771
0.848672914
0.612227091
0.140907759
0.93025033
0.55257952
0.006757016
0.295226718
0.192584359
0.217828162
1.340434943
0.310834204
2.739862161
2.591963755


#miles days* (lM)

-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
-9
## -------------------------------------------------------------------------- ##
## MARKER FOR END OF DATA FILE (eof1)                                          ##
## -------------------------------------------------------------------------- ##
42


