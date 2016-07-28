## —————————————————————————————————————————————————————————————————————————— ##
##                  DESIGN MATRIX FOR PARAMETER CONTROLS                      ##
##  Prior descriptions   Parameter values                                     ##                                                       ##
##  -0 uniform           (0,0)                                                ##
##  -1 normal            (p1=mu,p2=sig)                                       ##
##  -2 lognormal         (p1=log(mu),p2=sig)                                  ##
##  -3 beta              (p1=alpha,p2=beta)                                   ##
##  -4 gamma             (p1=alpha,p2=beta)                                   ##
## —————————————————————————————————————————————————————————————————————————— ##
## init   lower   upper    est  prior
## valu   bound   bound    phz   type     p1    p2   # PARAMETER              ##
## —————————————————————————————————————————————————————————————————————————— ##
  -1.05   -6.79    1.00      2      1    0.2  0.05   # log_natural_mortality
   7.50   -6.00   12.00      1      0      0     0   # log_rinit
   7.50   -6.00   12.00      1      0      0     0   # log_rbar
   7.50   -6.00   12.00      2      0      0     0   # log_ro
   1.50    0.00   12.00      2      0      0     0   # log_reck
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##                CONTROLS FOR TIME-VARYING MATURITY                          ##
## —————————————————————————————————————————————————————————————————————————— ##
## Phase for estimation if nMatBlocks > 1
  -2
## nMatBlocks
  1
## Year   
   2015
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##            CONTROLS FOR TIME-VARYING LN(NATURAL MORTALITY DEVS)            ##
## —————————————————————————————————————————————————————————————————————————— ##
## Phase for estimation if nMortBlocks > 1
  -2
## nMortBlocks
  2
## The terminal year of each block
  1999  2015
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##                    CONTROLS FOR SELECTIVITY PARAMETERS                     ##
## —————————————————————————————————————————————————————————————————————————— ##
## - Each selectivity block can have different functional forms based on selType
## - LEGEND:
##   - SelType = 1, logistic selectivity, 2 parameters.
##  nSelexblocks
    2
## —————————————————————————————————————————————————————————————————————————— ##
##  Gear  Sel     sel   sel   age   year  phz    | start end
##  Index Type    mu    sd    nodes nodes mirror | block block
## —————————————————————————————————————————————————————————————————————————— ##
    1     1       3.0   0.5   0     0      2        1980  2000
    1     1       5.0   0.3   0     0      2        2001  2015
## —————————————————————————————————————————————————————————————————————————— ##


## —————————————————————————————————————————————————————————————————————————— ##
##                        OTHER MISCELLANEOUS CONTROLS                        ##
## —————————————————————————————————————————————————————————————————————————— ##
## number of controls to read in.
2
## Value    # - Description
0.90718     # - Catch Scaler (convert from short tons to metric tons)
1           # - Condition on Catch = 0, Condition of Effort = 1

## EOF
999