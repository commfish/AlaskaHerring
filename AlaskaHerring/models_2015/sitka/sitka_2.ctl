## —————————————————————————————————————————————————————————————————————————— ##
##                  DESIGN MATRIX FOR PARAMETER CONTROLS                      ##
##  Prior descriptions   Parameter values                                     ##                                                                            ##
##  -0 uniform           (0,0)                                                ##
##  -1 normal            (p1=mu,p2=sig)                                       ##
##  -2 lognormal         (p1=log(mu),p2=sig)                                  ##
##  -3 beta              (p1=alpha,p2=beta)                                   ##
##  -4 gamma             (p1=alpha,p2=beta)                                   ##
## —————————————————————————————————————————————————————————————————————————— ##
##  init  lower   upper    est  prior
## value  bound   bound    phz   type     p1    p2   # PARAMETER              ##
## —————————————————————————————————————————————————————————————————————————— ##
   -1.05   -6.79    1.00      2      0  -1.05  0.05   # log_natural_mortality
    4.60   -6.00   12.00      1      0      0     0   # log_rinit
    5.60   -6.00   12.00      1      0      0     0   # log_rbar
    6.00   -6.00   12.00      2      0      0     0   # log_ro
    2.40    0.00   12.00      2      2   2.73  0.50   # log_reck
    6.25    0.00  200.00     -2      4   1.05  1.05   # precision = 1/(sigma_r^2)
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##                CONTROLS FOR TIME-VARYING MATURITY                          ##
## —————————————————————————————————————————————————————————————————————————— ##
## nMatBlocks
  1
## a50    a95     phz   terminalBlockYear
   4.5    7.0      -2           2015
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
  -2
## nMortBlocks, or Nodes in the case of cubic spline interpolation
  4
## The terminal year of each block
 1977 1998 2002 2015
## —————————————————————————————————————————————————————————————————————————— ##




## —————————————————————————————————————————————————————————————————————————— ##
##                    CONTROLS FOR SELECTIVITY PARAMETERS                     ##
## —————————————————————————————————————————————————————————————————————————— ##
## - Each selectivity block can have different functional forms based on selType
## - LEGEND:
##   - SelType = 1, logistic selectivity, 2 parameters.
##   - SelType = 2, logistic with 50% & 95% parameters
##  nSelexblocks
    3
## —————————————————————————————————————————————————————————————————————————— ##
##  Gear  Sel     sel   sel   age   year  phz    | start end
##  Index Type    mu    sd    nodes nodes mirror | block block
## —————————————————————————————————————————————————————————————————————————— ##
    1     1       3.0   0.5   0     0     -1        1971  1980
    1     1       5.0   0.3   0     0      2        1981  2000
    1     1       5.0   0.3   0     0      2        2001  2015
## —————————————————————————————————————————————————————————————————————————— ##


## —————————————————————————————————————————————————————————————————————————— ##
##                        OTHER MISCELLANEOUS CONTROLS                        ##
## —————————————————————————————————————————————————————————————————————————— ##
## number of controls to read in.
6
## Value    # # - Description
0.90718     # 1 - Catch Scaler (convert from short tons to metric tons)
1           # 2 - Condition on Catch = 0, Condition of Ft = 1
25000       # 3 - harvest threshold
0.10        # 4 - target harvest rate
20000       # 5 - threshold denominator
0.001       # 6 - standard deviation in natural mortality devs
## EOF
999