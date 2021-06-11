## Alaska Herring Assessment Model
Steve Martell.
May 12, 2016
Jun 29, 2016
Jul 20, 2016
Jul 29, 2016
----

## Introduction
This Repository documents changes made to the Age-structured Model for Alaska herring stocks, VERSION 0.1, Jan 2015.  The authors of the assessment are Sherri Dressel, Sara Miller, and Kray Van Kirk.  The code in VERSION 0.1 was developed by Peter Hulson <pete.hulson@noaa.gov> .

## Installation
You can obtain the source code via cloning this project using Git, or you may also download a zip file with the latest code. 


Go to the [Github code repository for Ham](https://github.com/seastateinc/AlaskaHerring), and click on the Clone or download button.

![Select the Download Zip option.](https://github.com/seastateinc/AlaskaHerring/blob/develop/docs/CloneZip.png)

I would recommend using Git and Github, as they serve as valuable tools for version control, and tracking code changes over time.  Its perfectly fine to just download the zip file; however, using Git is far more efficient for incorporating code changes and working in groups.

## R-scripts
The R-scripts directory has two important files:
	1. modelList.txt
	2. R4HAM.R

The modelList.txt file is a list of models that R4HAM.R will read into R.


----

## To Do List

- [x] Improve Numerical Stability:
	- [x] Parameter transformations to log-space
	- [x] Design matrix for estimated parameters (ival, lb, ub, phz, bayes)
	- [x] Develop option to condition the model on F and fit to catch.

- [ ] Code reorganization:
	- [x] HEAD::initializeModelParameters
	- [x] BIOL::initializeAgeSchedule
		- [x] BIOL::getNaturalMortality
		- [x] OBSM::getSelex
	- [x] BIOL::initializeStateVariables
	- [x] BIOL::updateStateVariables
	- [x] BIOL::calcSpawningStockRecruitment
	- [x] OBSM::calcAgeCompResiduals
	- [x] OBSM::calcEggMiledaySurveyResiduals
	- [x] OBSM::catchObservationModel
	- [x] STAT::calcObjectiveFunction
		- [x] STAT::penaltyFunctions
		- [x] STAT::negativeLogLikelihoods
		- [x] STAT::constraintFunctions
		- [x] STAT::calculateDIC
	- [ ] FORE::runForecast
		- [ ] FORE::ghlCalc
		- [ ] FORE::calcTAC

- [x] Control file reorganization
	- [x] add design matrix to control parameter bounds and phases. 
	- [x] add contrls for time varying maturity
	- [x] add contrls for time varying natural mortality rate deviations.
	- [x] add design matrix for selectivity parameter controls.
	- [x] add Miscellaneous controls for appending.


- [ ] Simulation model:
	- [x] Add command line argument to turn on simulation model.
	- [x] Add FUNCTION runSimulationModel
		- [x] Get and use True Parameter values
		- [x] Generate random variables | random number seed.
		- [ ] Population dynamics conditioned on catch & S-R curve.
		- [x] Cache observation model results into data variables.
		- [x] Add observation errors.
		- [x] Allow Alaska Herring Assessment to continue with parameter estimation.

## Reproducibility ----

Last updated: June 2020

```
 setting  value                       
 version  R version 3.5.3 (2019-03-11)
 os       Windows >= 8 x64            
 system   x86_64, mingw32             
 ui       RStudio                     
 language (EN)                        
 collate  English_United States.1252  
 ctype    English_United States.1252  
 tz       America/Anchorage           
 date     2020-06-11                  

- Packages -------------------------------------------------------------------------------------------
 package     * version  date       lib source                                 
 assertthat    0.2.1    2019-03-21 [1] CRAN (R 3.5.3)                         
 backports     1.1.5    2019-10-02 [1] CRAN (R 3.5.3)                         
 broom       * 0.5.2    2019-04-07 [1] CRAN (R 3.5.3)                         
 callr         3.3.2    2019-09-22 [1] CRAN (R 3.5.3)                         
 captioner   * 2.2.3    2015-07-16 [1] CRAN (R 3.5.3)                         
 cellranger    1.1.0    2016-07-27 [1] CRAN (R 3.5.3)                         
 cli           1.1.0    2019-03-19 [1] CRAN (R 3.5.3)                         
 coda        * 0.19-3   2019-07-05 [1] CRAN (R 3.5.3)                         
 colorspace    1.4-1    2019-03-18 [1] CRAN (R 3.5.3)                         
 cowplot     * 1.0.0    2019-07-11 [1] CRAN (R 3.5.3)                         
 crayon        1.3.4    2017-09-16 [1] CRAN (R 3.5.3)                         
 crosstalk     1.0.0    2016-12-21 [1] CRAN (R 3.5.3)                         
 data.table  * 1.12.2   2019-04-07 [1] CRAN (R 3.5.3)                         
 DBI           1.0.0    2018-05-02 [1] CRAN (R 3.5.3)                         
 dbplyr        1.4.2    2019-06-17 [1] CRAN (R 3.5.3)                         
 desc          1.2.0    2018-05-01 [1] CRAN (R 3.5.3)                         
 devtools      2.1.0    2019-07-06 [1] CRAN (R 3.5.3)                         
 digest        0.6.21   2019-09-20 [1] CRAN (R 3.5.3)                         
 dplyr       * 0.8.3    2019-07-04 [1] CRAN (R 3.5.3)                         
 ellipsis      0.3.0    2019-09-20 [1] CRAN (R 3.5.3)                         
 fansi         0.4.0    2018-10-05 [1] CRAN (R 3.5.3)                         
 FishLife    * 1.0.2    2019-08-25 [1] Github (james-thorson/FishLife@a63e0c8)
 forcats     * 0.4.0    2019-02-17 [1] CRAN (R 3.5.3)                         
 fs            1.3.1    2019-05-06 [1] CRAN (R 3.5.3)                         
 generics      0.0.2    2018-11-29 [1] CRAN (R 3.5.3)                         
 ggdendro      0.1-20   2016-04-27 [1] CRAN (R 3.5.3)                         
 ggformula   * 0.9.1    2019-01-12 [1] CRAN (R 3.5.3)                         
 ggplot2     * 3.2.1    2019-08-10 [1] CRAN (R 3.5.3)                         
 ggrepel     * 0.8.1    2019-05-07 [1] CRAN (R 3.5.3)                         
 ggstance    * 0.3.3    2019-08-19 [1] CRAN (R 3.5.3)                         
 ggthemes    * 4.2.0    2019-05-13 [1] CRAN (R 3.5.3)                         
 gh            1.0.1    2017-07-16 [1] CRAN (R 3.5.3)                         
 glue          1.3.1    2019-03-12 [1] CRAN (R 3.5.3)                         
 gridExtra   * 2.3      2017-09-09 [1] CRAN (R 3.5.3)                         
 gtable        0.3.0    2019-03-25 [1] CRAN (R 3.5.3)                         
 haven         2.2.0    2019-11-08 [1] CRAN (R 3.5.3)                         
 hms           0.5.3    2020-01-08 [1] CRAN (R 3.5.3)                         
 htmltools     0.3.6    2017-04-28 [1] CRAN (R 3.5.3)                         
 htmlwidgets   1.3      2018-09-30 [1] CRAN (R 3.5.3)                         
 httpuv        1.5.1    2019-04-05 [1] CRAN (R 3.5.3)                         
 httr          1.4.1    2019-08-05 [1] CRAN (R 3.5.3)                         
 jsonlite      1.6      2018-12-07 [1] CRAN (R 3.5.3)                         
 knitr       * 1.24     2019-08-08 [1] CRAN (R 3.5.3)                         
 labeling      0.3      2014-08-23 [1] CRAN (R 3.5.2)                         
 later         0.8.0    2019-02-11 [1] CRAN (R 3.5.3)                         
 lattice     * 0.20-38  2018-11-04 [2] CRAN (R 3.5.3)                         
 lazyeval      0.2.2    2019-03-15 [1] CRAN (R 3.5.3)                         
 leaflet       2.0.2    2018-08-27 [1] CRAN (R 3.5.3)                         
 lifecycle     0.1.0    2019-08-01 [1] CRAN (R 3.5.3)                         
 lubridate   * 1.7.4    2018-04-11 [1] CRAN (R 3.5.3)                         
 magrittr      1.5      2014-11-22 [1] CRAN (R 3.5.3)                         
 MASS          7.3-51.4 2019-04-26 [1] CRAN (R 3.5.3)                         
 Matrix      * 1.2-17   2019-03-22 [1] CRAN (R 3.5.3)                         
 memoise       1.1.0    2017-04-21 [1] CRAN (R 3.5.3)                         
 mime          0.7      2019-06-11 [1] CRAN (R 3.5.3)                         
 modelr        0.1.5    2019-08-08 [1] CRAN (R 3.5.3)                         
 mosaic      * 1.5.0    2019-01-12 [1] CRAN (R 3.5.3)                         
 mosaicCore    0.6.0    2018-06-24 [1] CRAN (R 3.5.3)                         
 mosaicData  * 0.17.0   2018-06-23 [1] CRAN (R 3.5.3)                         
 munsell       0.5.0    2018-06-12 [1] CRAN (R 3.5.3)                         
 nlme          3.1-137  2018-04-07 [2] CRAN (R 3.5.3)                         
 pillar        1.4.2    2019-06-29 [1] CRAN (R 3.5.3)                         
 pkgbuild      1.0.6    2019-10-09 [1] CRAN (R 3.5.3)                         
 pkgconfig     2.0.3    2019-09-22 [1] CRAN (R 3.5.3)                         
 pkgload       1.0.2    2018-10-29 [1] CRAN (R 3.5.3)                         
 plyr          1.8.4    2016-06-08 [1] CRAN (R 3.5.3)                         
 prettyunits   1.0.2    2015-07-13 [1] CRAN (R 3.5.3)                         
 processx      3.4.1    2019-07-18 [1] CRAN (R 3.5.3)                         
 promises      1.0.1    2018-04-13 [1] CRAN (R 3.5.3)                         
 ps            1.3.0    2018-12-21 [1] CRAN (R 3.5.3)                         
 purrr       * 0.3.4    2020-04-17 [1] CRAN (R 3.5.3)                         
 R2admb      * 0.7.16   2017-10-30 [1] CRAN (R 3.5.3)                         
 R6            2.4.0    2019-02-14 [1] CRAN (R 3.5.3)                         
 Rcpp          1.0.2    2019-07-25 [1] CRAN (R 3.5.3)                         
 readr       * 1.3.1    2018-12-21 [1] CRAN (R 3.5.3)                         
 readxl        1.3.1    2019-03-13 [1] CRAN (R 3.5.3)                         
 remotes       2.1.0    2019-06-24 [1] CRAN (R 3.5.3)                         
 reprex        0.3.0    2019-05-16 [1] CRAN (R 3.5.3)                         
 reshape2      1.4.3    2017-12-11 [1] CRAN (R 3.5.3)                         
 rfishbase     3.0.4    2019-06-27 [1] CRAN (R 3.5.3)                         
 rlang         0.4.5    2020-03-01 [1] CRAN (R 3.5.3)                         
 rprojroot     1.3-2    2018-01-03 [1] CRAN (R 3.5.3)                         
 rstudioapi    0.10     2019-03-19 [1] CRAN (R 3.5.3)                         
 rvest         0.3.5    2019-11-08 [1] CRAN (R 3.5.3)                         
 scales      * 1.0.0    2018-08-09 [1] CRAN (R 3.5.3)                         
 sessioninfo   1.1.1    2018-11-05 [1] CRAN (R 3.5.3)                         
 shape         1.4.4    2018-02-07 [1] CRAN (R 3.5.2)                         
 shiny         1.3.2    2019-04-22 [1] CRAN (R 3.5.3)                         
 stringi       1.4.3    2019-03-12 [1] CRAN (R 3.5.3)                         
 stringr     * 1.4.0    2019-02-10 [1] CRAN (R 3.5.3)                         
 testthat      2.2.1    2019-07-25 [1] CRAN (R 3.5.3)                         
 tibble      * 2.1.3    2019-06-06 [1] CRAN (R 3.5.3)                         
 tidyr       * 1.0.2    2020-01-24 [1] CRAN (R 3.5.3)                         
 tidyselect    0.2.5    2018-10-11 [1] CRAN (R 3.5.3)                         
 tidyverse   * 1.3.0    2019-11-21 [1] CRAN (R 3.5.3)                         
 usethis       1.5.1    2019-07-04 [1] CRAN (R 3.5.3)                         
 utf8          1.1.4    2018-05-24 [1] CRAN (R 3.5.3)                         
 vctrs         0.2.4    2020-03-10 [1] CRAN (R 3.5.3)                         
 withr         2.1.2    2018-03-15 [1] CRAN (R 3.5.3)                         
 xfun          0.9      2019-08-21 [1] CRAN (R 3.5.3)                         
 xml2          1.2.2    2019-08-09 [1] CRAN (R 3.5.3)                         
 xtable        1.8-4    2019-04-21 [1] CRAN (R 3.5.3)                         
 yaml          2.2.0    2018-07-25 [1] CRAN (R 3.5.3)                         

[1] C:/Users/jysullivan/Documents/R/win-library/3.5
[2] C:/Users/jysullivan/Documents/R/R-3.5.3/library
```
