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
  -1.05   -6.79    1.00      2      1   -0.7  0.05   # log_natural_mortality
   7.50   -6.00   12.00      1      0      0     0   # log_rinit
   7.50   -6.00   12.00      1      0      0     0   # log_rbar
   7.50   -6.00   12.00      2      0      0     0   # log_ro
   1.50    0.00   12.00      2      0      0     0   # log_reck
  -3.40    0.00    1.00     -2      0      0     0   # sigma_r
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##                CONTROLS FOR TIME-VARYING MATURITY                          ##
## —————————————————————————————————————————————————————————————————————————— ##
## nMatBlocks
  3
## a50    a95     phz   terminalBlockYear
   4.5    7.0       2           1999
   4.5    7.0       2           2003
   4.5    7.0       2           2015
## —————————————————————————————————————————————————————————————————————————— ##





## —————————————————————————————————————————————————————————————————————————— ##
##            CONTROLS FOR TIME-VARYING LN(NATURAL MORTALITY DEVS)            ##
## KEY:
##  Type: 1 = constant M
##  Type: 2 = interpolated using cubic spline.
## —————————————————————————————————————————————————————————————————————————— ##
## Type
  1
## Phase for estimation if nMortBlocks > 1
  2
## nMortBlocks
  3
## The terminal year of each block
  1999  2008 2015
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##                    CONTROLS FOR SELECTIVITY PARAMETERS                     ##
## —————————————————————————————————————————————————————————————————————————— ##
## - Each selectivity block can have different functional forms based on selType
## - LEGEND:
##   - SelType = 1, logistic selectivity, 2 parameters.
##  nSelexblocks
    1
## —————————————————————————————————————————————————————————————————————————— ##
##  Gear  Sel     sel   sel   age   year  phz    | start end
##  Index Type    mu    sd    nodes nodes mirror | block block
## —————————————————————————————————————————————————————————————————————————— ##
    1     1       3.0   0.5   0     0      2        1988  2015
## —————————————————————————————————————————————————————————————————————————— ##


## —————————————————————————————————————————————————————————————————————————— ##
##                        OTHER MISCELLANEOUS CONTROLS                        ##
## —————————————————————————————————————————————————————————————————————————— ##
## number of controls to read in.
5
## Value    # - Description
0.90718     # - Catch Scaler (convert from short tons to metric tons)
1           # - Condition on Catch = 0, Condition of Ft = 1
5000        # - harvest threshold
0.20        # - target harvest rate
5000        # - threshold denominator
## EOF
999