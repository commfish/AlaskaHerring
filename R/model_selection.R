
# Model selection for HER
# Jane Sullivan (jane.sullivan1@alaska.gov)
# Updated January 2019

# Pseudocode
# 1. User inputs the year breaks
# 2. Generate time block combos and permutations across the 3 time varying
# parameters
# 3. Remove unwanted models that violate user rules (e.g. selectivity can only
# change within the same years or less than maturity)
# 4. Run models in loop
# 5. Create a summary of models that did not converge
# 6. In final loop create a subdirectory of figures and extract AIC and summary
# statistics for model selection.

# Data/Libraries ----

library(tidyverse)  # data processing
library(gtools)     # for the permutations function
library(R2admb)     # running ADMB through R
library(data.table)
source("R/helper.R")
source("R/figure_fxns.R")

# User inputs ----

# Start year, last year, and breaks
syr <- 1980
lyr <- 2017
breaks <- c(1998, 2014) # last year of time block, defined by PDO or other method

# 1. Time block combos ----

# Reduce all combinations into possible ones, define each combo by the last year
# of each time block
allcombos <- do.call(c, lapply(seq_along(c(breaks, lyr)), combn, x = c(breaks, lyr), simplify = FALSE))
combos <- lapply(allcombos, function(x) max(x) == lyr)
combos <- allcombos[unlist(combos)]

# Make each combination a df with a column for start and end year of each time block 
blocks <- lapply(combos, function(x) data.frame(lyr = x))
blocks <- lapply(blocks, function(x) {x$syr = c(syr, head(x$lyr, -1) + 1); return(x)})
blocks <- lapply(blocks, function(x) x[c('syr', 'lyr')])

# Number of time blocks in each combo
n_blocks <- lapply(blocks, nrow)

# 2. Permuations by parameter ----

# Permutations of the time block combinations over the 3 time-varying parameters
# (always in the following order: survival, maturity, selectivity) 
perms <- gtools::permutations(n = length(combos), # size of source vector 
                              r = 3, # number of time-varying parameters
                              v = 1:length(combos), # source vector (each number represents a combo)
                              repeats.allowed = TRUE)
perms <- as.data.frame(perms)
names(perms) <- c("survival", "maturity", "selectivity")

# Function that applies rule that selectivity can only change within the same
# years or less than maturity). This eliminates all the extra models, which
# makes the total run time much faster. If you don't want to apply this rule,
# don't run the function
selectivity_tst <- function(perms = perms,
                            # If debug = TRUE the tst list will print out. see
                            # debug=TRUE documentation below for details
                            debug = FALSE) {
  
  tst <- list()
  for(i in 1:length(combos)) {
    vec <- vector() # temporary vector
    for(j in 1:length(combos)) {
      # If the jth combo is contained within the ith combo it passes the test
      if(all(combos[[j]] %in% combos[[i]])) { vec[j] <- j }
    }
    tst[[i]] <- vec[!is.na(vec)]
  }
  
  # debug = TRUE prints out the tst list:
  # Each element of this list is a maturity combo, and the vector of combos are the
  # acceptable selectivity combos Examples: If maturity is the first combo,
  # selectivity can only be combo 1. If maturity is combo 3, selectivity can by
  # combo 1 or 3
  if(debug) {print(tst)}

  # Filter these out of our permutation df
  out <- list()
  for(i in 1:length(combos)) {
    out[[i]] <- perms %>% filter(maturity == i & selectivity %in% tst[[i]])
  }
  perms <- do.call(rbind, out) # collapses list into df
  return(perms)
}

# Remove unwanted model runs where selectivity blocks occur outside of maturity
# blocks
perms <- selectivity_tst(perms)
# perms <- selectivity_tst(perms, debug = TRUE) # debug mode

# 3. Write ctl ----

# Function writes maturity block section in ctl file
write_maturity_blocks <- function(i) {
  combo <- perms[i,'maturity'] # Extract relevant combination
  out <- list()
  for(i in 1:n_blocks[[combo]]) {
    tmp <- blocks[[combo]]$lyr[i]
    # user-defined starting values for a50, a95 and estimation phase
    out[[i]] <- paste0("4.5       7.0       2       ", tmp, collapse = "     ")
  }
  txt <- paste(unlist(out), sep = "\n")
  return(txt)
}

# Function writes survival block section in ctl file
write_survival_blocks <- function(i) {
  combo <- perms[i,'survival'] # Extract relevant combination
  txt <- paste(blocks[[combo]]$lyr, collapse = " ")
  return(txt)
}

# Function writes selectivity block section in ctl file
write_selectivity_blocks <- function(i) {
  combo <- perms[i,'selectivity'] # Extract relevant combination
  out <- list()
  for(i in 1:n_blocks[[combo]]) {
    syr_tmp <- blocks[[combo]]$syr[i]
    lyr_tmp <- blocks[[combo]]$lyr[i]
    # user-defined starting values for gear index, sel type, sel mu and sd, age
    # nodes, year nodes, and the phase for estimation
    out[[i]] <- paste("1     1     3.0     0.5     0     0     2     ", 
                       syr_tmp, lyr_tmp, sep = "     ", collapse = "    ")
  }
  txt <- paste(unlist(out), sep = "\n")
  return(txt)
}

# Function that writes .ctl file
write_ctl <- function(i) {
ctl <- c(
"
## —————————————————————————————————————————————————————————————————————————— ##
##                  DESIGN MATRIX FOR PARAMETER CONTROLS                      ##
##  Prior descriptions   Parameter values                                     ##                                                                            ##
##  -0 uniform           (0,0)                                                ##
##  -1 normal            (p1=mu,p2=sig)                                       ##
##  -2 lognormal         (p1=log(mu),p2=sig)                                  ##
##  -3 beta              (p1=alpha,p2=beta)                                   ##
##  -4 gamma             (p1=alpha,p2=beta)                                   ##
## —————————————————————————————————————————————————————————————————————————— ##
##  init  lower   upper    est  prior                                         ##
## value  bound   bound    phz   type     p1    p2   # PARAMETER              ##
## —————————————————————————————————————————————————————————————————————————— ##
-1.05   -6.79    1.00      2      0  -1.05  0.05   # log_natural_mortality  (original, p1 too high?)
## -0.7985  -5.00    5.00     1      1  -0.7985  0.4  # log_natural_mortality (same as iscam)
4.60   -6.00   12.00      1      0      0     0   # log_rinit
5.60   -6.00   12.00      1      0      0     0   # log_rbar
6.00   -6.00   12.00      2      0      0     0   # log_ro
2.40    0.00   12.00      2      2   2.73  0.50   # log_reck
4.00    0.00  200.00     -2      4   1.05  1.05   # precision = 1/(sigma_r^2), SM=6.25
## 50.25    0.00  200.00    -2      4   1.05  1.05   # precision = 1/(sigma_r^2) (for simulations with no process error)

## —————————————————————————————————————————————————————————————————————————— ##

## —————————————————————————————————————————————————————————————————————————— ##
##                CONTROLS FOR TIME-VARYING MATURITY                          ##
## —————————————————————————————————————————————————————————————————————————— ##
## nMatBlocks",
n_blocks[[perms[i,'maturity']]],
"## a50    a95     phz   terminalBlockYear",
write_maturity_blocks(i),
"## —————————————————————————————————————————————————————————————————————————— ##

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
## nMortBlocks, or Nodes in the case of cubic spline interpolation",
n_blocks[[perms[i,'survival']]],
"## The terminal year of each block",
write_survival_blocks(i),
"## ————————————————————————————————————————————————————————————————————————— ##

## —————————————————————————————————————————————————————————————————————————— ##
##                    CONTROLS FOR SELECTIVITY PARAMETERS                     ##
## —————————————————————————————————————————————————————————————————————————— ##
## - Each selectivity block can have different functional forms based on selType
## - LEGEND:
##   - SelType = 1, logistic selectivity, 2 parameters.
##   - SelType = 2, logistic with 50% & 95% parameters
##  nSelexblocks",
n_blocks[[perms[i,'selectivity']]],
"## —————————————————————————————————————————————————————————————————————————— ##
##  Gear  Sel     sel   sel   age   year  phz    | start end
##  Index Type    mu    sd    nodes nodes mirror | block block
## —————————————————————————————————————————————————————————————————————————— ##",
write_selectivity_blocks(i),
"## ——————————————————————————————————————————————————————————————————————————— ##

## ——————————————————————————————————————————————————————————————————————————— ##
##                        OTHER MISCELLANEOUS CONTROLS                         ##
## ——————————————————————————————————————————————————————————————————————————— ##
## number of controls to read in.
8
## Value    # # - Description
0.90718     # 1 - Catch Scaler (convert from short tons to metric tons)
0           # 2 - Condition on Catch = 0, Condition of Ft = 1
25000       # 3 - harvest threshold
0.2         # 4 - target harvest rate
20000       # 5 - threshold denominator
0.075       # 6 - standard deviation in natural mortality devs SM=0.001
1.00        # 7 - sd in recruitment deviations in all phases of estimate up until the last
5.00        # 8 - sd in recruitment deviation in the final phase of estimation
## EOF
999
")
# Write it to file
write.table(ctl, file = "sitka.ctl", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# 4. Directory setup ----
root_dir <- getwd() # project root
tpl_dir <- file.path(root_dir, "HER") # location of tpl
dir.create("Model_selection") # new folder for model selection
modsel_dir <- file.path(root_dir, "Model_selection")
files <- list.files(file.path(tpl_dir)) # searchable file list
tpl <- file.path(tpl_dir, files[which(grepl(".tpl", files))]) # find .tpl
name <- strsplit(files[grepl("tpl",files)],".tpl")[[1]][1] # program name
dat <- file.path(tpl_dir, files[which(grepl(".dat", files))]) # fine .dat files (two)
std <- paste0(name,".std") # used to check model convergence
exe <- paste0(name,".exe") # executible file so it doesn't need to be recompiled each time

# 5. Run models ----

# Function to format the time blocks for each parameter, regardless of how many
# blocks there are
write_blks <- function(i, param = 'survival') {
  # extract related start and end year for the time block combo
  blks <- blocks[[perms[i, param]]] 
  out <- list()  # store output
  for(j in 1:nrow(blks)) {
    out[[j]] <- paste0(blks$syr[j], "-", blks$lyr[j]) } # formating
  txt <- paste(unlist(out), collapse = "; ") 
  return(txt)
}

# Convergence check: 0 = not converged, 1 = converged
check_converge <- matrix(ncol = 7, nrow = nrow(perms)) 

for(i in 1:nrow(perms)){
  
  # Model naming convention: model combination # followed by number of time
  # blocks for natural survival, maturity, and selectivity
  model_name <- paste0("model.", i, "_",
                       n_blocks[[perms[i,'survival']]],
                       n_blocks[[perms[i,'maturity']]],
                       n_blocks[[perms[i,'selectivity']]])
  
  # Create model-specific directory
  model_dir <- file.path(modsel_dir, model_name)
  dir.create(model_dir, showWarnings = FALSE)
  setwd(model_dir)
  
  # Create new sitka.ctl file, copy tpl, exe, and dat files to model-specific
  # directory
  write_ctl(i) 
  file.copy(from = tpl, to = model_dir, overwrite = TRUE)	
  file.copy(from = file.path(tpl_dir, exe), to = model_dir, overwrite = TRUE)	
  file.copy(from = dat[1], to = model_dir, overwrite = TRUE)	
  file.copy(from = dat[2], to = model_dir, overwrite = TRUE)	
  
  # remove pre-existing std file if it exists and run admb
  if (file.exists(std)) {file.remove(std, showWarnings = FALSE)}
  run_admb(name, verbose = TRUE)
  
  # store formatted convergence results
  check_converge[i,1] <- i
  check_converge[i,2] <- model_name
  check_converge[i,3] <- paste0("Surv(", n_blocks[[perms[i,'survival']]],")", 
                                " Mat(", n_blocks[[perms[i,'maturity']]],")", 
                                " Sel(", n_blocks[[perms[i,'selectivity']]],")")
  check_converge[i,4] <- write_blks(i, param = 'survival')
  check_converge[i,5] <- write_blks(i, param = 'maturity')
  check_converge[i,6] <- write_blks(i, param = 'selectivity')
  # check for convergence (1 = converged, 0 = not converged)
  if (file.exists(std)) {check_converge[i,7] <- 1} 
  else {check_converge[i,7] <- 0 }
}

# 6. Check convergence ----
check_converge <- as.data.frame(check_converge)
names(check_converge) <- c("Model", "Folder_name", "Structure", "Survival", 
                           "Maturity", "Selectivity", "Convergence")

# Save relevant files for later use 
setwd(modsel_dir)

if(exists("check_converge")) {
  write_csv(check_converge, "convergence.csv")
} else {
  check_converge <- read_csv("convergence.csv") }

# Models that did not converge - as of 2019-01-31 only the models with a single
# time block for maturity did not converge
filter(check_converge, Convergence == 0)

# Filter out models that did not converge
filter(check_converge, Convergence == 1) -> selection_summary

# Save output for final model selection summaries
selection_summary %>% 
  mutate(AIC = NA,
         max_grad = NA,
         GHL = NA,
         nopar = NA) -> selection_summary

# 7. Figures and report file ----

# total number of years in the data vs model and model start year index, used
# for indexing output
n_dat_yrs <- D[["dat_nyr"]] - D[["dat_syr"]] + 1
n_mod_yrs <- D[["mod_nyr"]] - D[["mod_syr"]] + 1
syr_index <- D[["mod_syr"]] - D[["dat_syr"]] + 1

# Function that formats matrices in the format year x age
format_matrix <- function(D = D, 
                          # Name of year x age matrix you want formatted
                          matrix = 'Nij', 
                          # 'years' includes forecast year, 'year' doesn't
                          yr = 'years') {
  Year <- data.frame(Year = D[[yr]])
  mat <- as.data.frame(D[[matrix]])
  names(mat) <- c(3:7, "8+")
  mat <- cbind(Year, mat)
  mat <- mat %>% filter(Year <= D[['mod_nyr']] & Year >= D[['mod_syr']])
  return(mat)
}

for(i in 1:length(selection_summary$Folder_name)){
  
  # Navigate to the correct directory
  fn <- selection_summary$Folder_name[i]
  model_dir <- file.path(modsel_dir, fn)
  setwd(model_dir)
  
  D <- read_admb(name) # Results
  P <- read_fit(name)  # Diagnostics
  
  # Create directory for figures and results
  results <- file.path(model_dir, "Results")
  dir.create(results, showWarnings = FALSE)
  
  # Store salient results for selection summary
  selection_summary[i,"AIC"] <- 2 * P[["nlogl"]] + 2 * P[["nopar"]]
  selection_summary[i,"max_grad"] <- P[["maxgrad"]]
  selection_summary[i, "GHL"] <- D[["ghl"]] # Already in short tons
  selection_summary[i, "nopar"] <- P[["nopar"]]
  
  # Generate figures
  plot_recruitment(D, results, showRicker = FALSE) # showRicker=TRUE Ricker estimates in top graph desired
  plot_eggdep(D, results, 
              # Show ootstrap confidence intervals
              bootstrap = TRUE,
              # Show 95% credibility interval from posterior samples (not
              # possible unless you run the MCMC first)
              credibility = FALSE, 
              # Cut off extreme value in 2008
              remove_outlier = TRUE)
  plot_survival(D, results)
  plot_matsel(D, results)
  # generates 4 figs (data bubble plot, residuals, & barplots with data and fit)
  plot_agecomps(D, results) 
  plot_biomass(D, results)
  
  # Write report
  
  # Forecast values
  forec <- data.frame(age = c(3:7, "8+"),
                      nj = D[['fore_nj']],
                      matnj = D[['fore_matnj']],
                      prop_matnj = D[['fore_matnj']] / sum(D[['fore_matnj']]),
                      waa = D[['data_sp_waa']][n_dat_yrs, 2:7],
                      surv = exp(-D[['Mij']][n_mod_yrs, ]),
                      mat = D[['mat']][n_mod_yrs, ],
                      sel = D[['Sij']][n_mod_yrs, ])
  
  names(forec) <- c('Age', 'Mature and immature numbers-at-age forecast (millions)',
                    'Mature numbers-at-age (millions)', 'Proportion mature numbers-at-age',
                    'Weight-at-age used in forecast (g)', 'Survival', 'Maturity', 'Selectivity')
  
  # Derived time series
  ts <- data.frame(Year = D[['year']],
                   matB = D[['mat_B']] / 0.90718, # convert to short tons
                   catch = D[['data_catch']][syr_index:n_dat_yrs, 2] / 0.90718, 
                   rec = D[['Nij']][1:n_mod_yrs, 1]) %>% 
    mutate(spB = matB - catch)  # Assumption that 100% catch is mature
  
  ts <- ts %>% 
    select(Year, `Mature biomass (tons)` = matB,
           `Spawning biomass (tons)` = spB, `Catch (tons)` = catch,
           `Mature and immature age-3 abundance (recruitment in millions)` = rec)
  
  # Format year x age matrices
  Nij <- format_matrix(D, 'Nij') # Total numbers-at-age 
  mat_Nij <- format_matrix(D, 'mat_Nij') # Mature numbers-at-age
  Cij <- format_matrix(D, 'Cij') # Catch in numbers-at-age
  # Spawning numbers-at-age with assumption that 100% of the catch is mature
  sp_Nij <- mat_Nij - Cij 
  sp_Nij <- sp_Nij %>% mutate(Year = D[['year']])
  # Estimated proportion of each age class that is caught (e.g. number of age-3
  # fish caught divided by the number of age-3 mature at age)
  Pij <- Cij / mat_Nij
  Pij <- Pij %>% mutate(Year = D[['year']])
  mat_Bij <- format_matrix(D, 'mat_Bij') # Mature biomass-at-age
  # Catch biomass-at-age
  catch_Bij <- D[['data_cm_waa']][syr_index:n_dat_yrs, ] * Cij
  catch_Bij <- catch_Bij %>% mutate(Year = D[['year']])
  # Spawning biomass-at-age with assumption at 100% of the catch is mature
  sp_Bij <- mat_Bij - catch_Bij
  # Surivival
  surv <- as.data.frame(exp(-D[['Mij']]))
  names(surv) <-  c(3:7, "8+")
  surv <- cbind(data.frame(Year = D[['year']]), surv)
  mat <- format_matrix(D, 'mat', yr = 'year') # Maturity
  sel <- format_matrix(D, 'Sij') # Selectivity
  # Observed spawning age comps
  sp_comp <- as.data.frame(D[['data_sp_comp']]) 
  names(sp_comp) <- c('Year', 3:7, "8+")
  sp_comp <- sp_comp %>% filter(Year <= D[['mod_nyr']] & Year >= D[['mod_syr']])
  # Observed catch age comps
  catch_comp <- as.data.frame(D[['data_cm_comp']]) 
  names(catch_comp) <- c('Year', 3:7, "8+")
  catch_comp <- catch_comp %>% filter(Year <= D[['mod_nyr']] & Year >= D[['mod_syr']])
  
  rep <- c(paste0("MODEL RESULTS", "\n",
                  "\n",
                  "Report produced by model_selection.R", "\n",
                  "Contact: Sherri.Dressel@alaska.gov or Jane.Sullivan1@alaska.gov", "\n",
                  paste0(selection_summary$Folder_name[i]),  "\n",
                  "\n",
                  "Model structure", "\n",
                  "Survival time blocks:,", selection_summary$Survival[i], "\n",
                  "Maturity time blocks:,", selection_summary$Maturity[i], "\n",
                  "Selectivity time blocks:,", selection_summary$Selectivity[i], "\n",
                  "\n",
                  "Model diagnostics", "\n",
                  "Number of parameters:,", round(P[['nopar']], 0),  "\n",
                  "Negative log likelihood:,", round(P[['nlogl']], 1),  "\n",
                  "Maximum gradient compoent:,", P[['maxgrad']], "\n",
                  "\n",
                  "Forecast", "\n",
                  "Mature biomass forecast (tons):,", round(D[['fore_matb']] / 0.90718, 1),  "\n",
                  "Threshold placeholder, TODO",  "\n",
                  "GHL (tons):,", round(D[['ghl']], 1),  "\n",
                  "Harvest rate:, 0.2", "\n",
                  "\n",
                  "Forecast values by age")
  )
  
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n")
  write.table(forec, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Catch and biomass time series"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(ts, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Mature and immature numbers-at-age (millions)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(Nij, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Mature numbers-at-age (millions)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(mat_Nij, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Spawning numbers-at-age (millions)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(sp_Nij, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Catch in numbers-at-age (millions)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(Cij, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Proportion of each age class that is caught from mature population (e.g. number of age-3 fish caught divided by the number of age-3 mature at age)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(Pij, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Mature biomass-at-age (tons)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(mat_Bij, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Spawning biomass-at-age (tons)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(sp_Bij, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Catch in biomass-at-age (tons)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(catch_Bij, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Survival"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(surv, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Maturity"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(mat, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Selectivity"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(sel, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Observed spawning age composition (cast net)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(sp_comp, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  rep <- c(paste0("\n", "Observed commercial age composition (spring seine)"))
  write.table(rep, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
  write.table(catch_comp, file = paste0(results, "/Report.csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
  
  }

# 8. Final AIC summary ----
selection_summary %>% 
  mutate(deltaAIC = abs(min(selection_summary$AIC) - AIC)) %>% 
  select(Folder_name, Structure, deltaAIC, AIC, GHL, nopar, max_grad, Survival, Maturity, Selectivity) %>% 
  arrange(deltaAIC) -> selection_summary

write_csv(selection_summary, paste0(modsel_dir, "/AIC_selection_summary.csv"))
