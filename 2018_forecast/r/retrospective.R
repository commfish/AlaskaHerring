# Retrospective analysis
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-06-03

# This script runs a sensitivity analysis on the HER model, generates figures,
# and calculates Mohn's and Woods Hole Rho.

# Libraries/source files and user inputs ----

library(R2admb) # run ADMB in R
library(tidyverse)
library(data.table)

# Forecast year, model start and last yrs
YEAR <- 2018
syr <- 1980 
lyr <- YEAR - 1

# IMPORTANT: this will create a subdirectory and allow you to run a
# retrospective on multiple model versions for comparison
MODEL_VERSION <- "HER_bestLS"   # HER with the "best" LS parameterization

# For MODEL_VERSION, what are the time blocks? This information is needed
# to accomodate short time blocks. Vector is the last year in each time block
survival_blk <- c(1998, 2014, lyr) # Natural mortality/survival
maturity_blk <- c(2014, lyr) # Maturity
selectivity_blk <- c(lyr) # Selectivity

source(paste0(YEAR, "_forecast/r/helper.r"))

# Directory set up ----

root_dir <- getwd() # project root
tpl_dir <- file.path(root_dir, paste0(YEAR, "_forecast/admb/HER_bestLS")) # location of tpl
results_dir <- file.path(root_dir, paste0(YEAR, "_forecast/results")) # results
dir.create(paste0(results_dir, "/retrospective/", MODEL_VERSION)) # new folder for retrospective
run_dir <- file.path(results_dir, paste0("retrospective/", MODEL_VERSION)) # folder to run the retrospective in
files <- list.files(file.path(tpl_dir)) # searchable file list
tpl <- file.path(tpl_dir, files[which(grepl(".tpl", files))]) # find .tpl
ctl <- file.path(tpl_dir, files[which(grepl(".ctl", files))]) # find .ctl
name <- strsplit(files[grepl("tpl",files)],".tpl")[[1]][1] # program name
dat <- file.path(tpl_dir, files[which(grepl(".dat", files))]) # fine .dat files (two)
std <- paste0(name,".std") # used to check model convergence
exe <- paste0(name,".exe") # executible file so it doesn't need to be recompiled each time

# Adjust control file to accomodate short time blocks ----

# There is an issue accomodating short time blocks at the end of the model time
# series. Sherri Dressel decided 2019-05-21 to apply the following rule: for
# retrospective peels that results in time blocks <= 2 yrs, remove time block.
# See https://github.com/commfish/AlaskaHerring/issues/4 for documentation.

# The following functions adjust the time blocks based on this rule and rewrite
# the control file for each retrospective peel. This framework currently only 
# works for <= 3 time blocks and the user must check that the fixed values in
# the write_ctl() fxn match those in the MODEL_VERSION's ctl.

retro_rule <- function(x) {
  
  # WARNING: This function only works for <= 3 time blocks.
  if(length(x) > 3) {return(print("STOP! This function only works for <= 3 time blocks. You have 4 or more time blocks and will need to extend this function!"))}
  
  # If there is only 1 time block no changes are needed
  if(length(x) == 1) {return(x)} 
  
  # Retrospective rule: for peels that results in time blocks <= 2 yrs,
  # remove the time block.
  
  # For 2 time blocks:
  if(length(x) == 2) {
    last2blks <- tail(x, 2)
    if(last2blks[2] - last2blks[1] <= 2) {
      x <- c(x[2]) # removing last time block
      return(x)}  else {return(x)}
  }
  
  # For 3 time blocks:
  if(length(x) == 3) {
    last2blks <- tail(x, 2)
    if(last2blks[2] - last2blks[1] <= 2) {
      x <- c(x[1], x[3]) # removing last time block
      return(x)} else {return(x)}
  } 
}

adjust_blks <- function(x) {
  
  # Create new time block vectors based on peel
  new_lyr <- lyr - retro[i]
  x <- replace(x, length(x), new_lyr)
  
  # Adjust time blocks based on retrospective rule on short time blocks
  x <- retro_rule(x)
  
  # Replace new_lyr with original lyr
  x <- replace(x, x == new_lyr, lyr)
  
  # Reformat
  x <- data.frame(lyr = x)
  x$syr <- c(syr, head(x$lyr, -1) + 1)
  x <- x[c('syr', 'lyr')]
  
  return(x)
}

# Function writes maturity block section in ctl file
write_maturity_blocks <- function(i) {
  out <- list()
  for(i in 1:length(new_maturity_blk$lyr)) {
    tmp <- new_maturity_blk$lyr[i]
    # user-defined starting values for a50, a95 and estimation phase
    out[[i]] <- paste0("4.5       7.0       2       ", tmp, collapse = "     ")
  }
  txt <- paste(unlist(out), sep = "\n")
  return(txt)
}

# Function writes selectivity block section in ctl file
write_selectivity_blocks <- function(i) {
  out <- list()
  for(i in 1:length(new_selectivity_blk$lyr)) {
    syr_tmp <- new_selectivity_blk$syr[i]
    lyr_tmp <- new_selectivity_blk$lyr[i]
    # user-defined starting values for gear index, sel type, sel mu and sd, age
    # nodes, year nodes, and the phase for estimation
    out[[i]] <- paste("1     1     5.0     0.3     0     0     2     ", 
                      syr_tmp, lyr_tmp, sep = "     ", collapse = "    ")
  }
  txt <- paste(unlist(out), sep = "\n")
  return(txt)
}

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
    ##  init  lower   upper    est  prior
    ## value  bound   bound    phz   type     p1    p2   # PARAMETER              ##
    ## —————————————————————————————————————————————————————————————————————————— ##
    -1.05   -6.79    1.00      2      0  -1.05  0.05  # log_natural_mortality
    4.60   -6.00   12.00      1      0      0     0   # log_rinit
    5.60   -6.00   12.00      1      0      0     0   # log_rbar
    6.00   -6.00   12.00      2      0      0     0   # log_ro
    2.40    0.00   12.00      2      2   2.73  0.50   # log_reck
    4.00    0.00  200.00     -2      4   1.05  1.05   # precision = 1/(sigma_r^2), SM=6.25
    ##  50.25    0.00  200.00    -2      4   1.05  1.05   # precision = 1/(sigma_r^2)
    
    ## —————————————————————————————————————————————————————————————————————————— ##

    ## —————————————————————————————————————————————————————————————————————————— ##
    ##                CONTROLS FOR TIME-VARYING MATURITY                          ##
    ## —————————————————————————————————————————————————————————————————————————— ##
    ## nMatBlocks",
    length(new_maturity_blk$lyr),
    "
    ## a50    a95     phz   terminalBlockYear",
    write_maturity_blocks(i),
    "
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
    ## nMortBlocks, or Nodes in the case of cubic spline interpolation",
    length(new_survival_blk$lyr),
    "
    ## The terminal year of each block",
    paste(new_survival_blk$lyr, collapse = " "),
    "
    ## —————————————————————————————————————————————————————————————————————————— ##

    ## —————————————————————————————————————————————————————————————————————————— ##
    ##                    CONTROLS FOR SELECTIVITY PARAMETERS                     ##
    ## —————————————————————————————————————————————————————————————————————————— ##
    ## - Each selectivity block can have different functional forms based on selType
    ## - LEGEND:
    ##   - SelType = 1, logistic selectivity, 2 parameters.
    ##   - SelType = 2, logistic with 50% & 95% parameters
    ##  nSelexblocks",
    length(new_selectivity_blk$lyr),
    "
    ## —————————————————————————————————————————————————————————————————————————— ##
    ##  Gear  Sel     sel   sel   age   year  phz    | start end
    ##  Index Type    mu    sd    nodes nodes mirror | block block
    ## —————————————————————————————————————————————————————————————————————————— ##",
    write_selectivity_blocks(i),
    "
    ## —————————————————————————————————————————————————————————————————————————— ##
    
    ## —————————————————————————————————————————————————————————————————————————— ##
    ##                        OTHER MISCELLANEOUS CONTROLS                        ##
    ## —————————————————————————————————————————————————————————————————————————— ##
    ## number of controls to read in.
    8
    ## Value    # # - Description
    0.90718     # 1 - Catch Scaler (convert from short tons to metric tons)
    0           # 2 - Condition on Catch = 0, Condition of Ft = 1
    25000       # 3 - harvest threshold
    0.2         # 4 - target harvest rate
    20000       # 5 - threshold denominator
    0.08        # 6 - standard deviation in natural mortality devs SM=0.001
    1.00        # 7 - sd in recruitment deviations in all phases of estimate up until the last
    5.00        # 8 - sd in recruitment deviation in the final phase of estimation
    ## EOF
    999
    
    ")
  # Write it to file
  write.table(ctl, file = "sitka.ctl", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Run retrospective ----

retro <- 0:10 # number of peels
check_converge <- list() # convergence results: 1 = converged, 0 = not converged

# Write ctl with updated time blocks - NOTE: You will need to check that the
# fixed values in this ctl function are the same for the MODEL_VERSION
for(i in 1:length(retro)){
  
  retro_dir <- file.path(run_dir, paste0("retro_", retro[i]))
  dir.create(retro_dir, showWarnings = FALSE)
  setwd(retro_dir)
  
  # Create new time block vectors based on peel
  new_survival_blk <- adjust_blks(survival_blk)
  new_maturity_blk <- adjust_blks(maturity_blk)
  new_selectivity_blk <- adjust_blks(selectivity_blk)
  
  # Create new sitka.ctl file, copy tpl, exe, and dat files to peel-specific
  # directory
  write_ctl(i)
  file.copy(from = tpl, to = retro_dir, overwrite = TRUE)	
  file.copy(from = file.path(tpl_dir, exe), to = retro_dir, overwrite = TRUE)	
  file.copy(from = dat[1], to = retro_dir, overwrite = TRUE)	
  file.copy(from = dat[2], to = retro_dir, overwrite = TRUE)	
  
  # remove pre-existing std file if it exists
  if (file.exists(std)) {file.remove(std, showWarnings = FALSE)}
  
  # run admb, set seed at retro value
  run_admb(name, extra.args=paste("-retro", retro[i]))
  
  # check for minimal convergence - does admb produce std file? is max grad component in par file < 1e-3? (1 = converge, 0 = not converged)
  if (file.exists(std)) {
    check_converge[[i]] <- 1
  } else { check_converge[[i]] <- 0 }
  
  mgc <- readMGC(name)
  
  if (mgc < 1e-3) {
    check_converge[[i]] <- 1
  } else {check_converge[[i]] <- 0 } 
  
}

# Check to see which models converges and which didn't
check_converge

# Graphics ----

# Empty vectors
SSB_ls <- list()

# Read results (only for converged models)
converged <- data.frame(retro = retro, converge = unlist(check_converge)) %>% 
  filter(converge == 1)

for(i in 1:length(converged$retro)){
  
  j <- pull(slice(converged, i), retro) # slice pulls row i
  retro_dir <- file.path(run_dir, paste0("retro_", j)) 
  setwd(retro_dir)
  
  D <- read_admb(name)
  # D <- read_rep(name)
  
  SSB_ls[[i]] <- data.frame(Year = D[["year"]], 
                         spB = D[["sp_B"]] / 0.90718, # convert to short tons
                         retro = paste0("retro_", j)) 
}

# collapse list
SSB <- do.call(rbind, SSB_ls)

axisx <- tickr(SSB, Year, 5)
# Generic retrospective plot
ggplot(SSB, aes(x = Year, y = spB, colour = retro, group = retro)) +
  geom_line(size = 1) +
  scale_color_grey(guide = FALSE) +
  labs(x = NULL, y = "Spawning biomass (tons)\n") +
  scale_y_continuous(label = comma) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) -> retro1

# Per Clark et al 2012, also show relative percent difference between peel and
# terminal year estimates
SSB %>% 
  group_by(retro) %>% 
  full_join(SSB %>% 
              filter(retro == "retro_0") %>% 
              select(Year, term_spB = spB),
            by = "Year") %>% 
  mutate(perc_diff = ((term_spB - spB) / term_spB) * 100) -> SSB

ggplot(data = SSB) +
  geom_hline(yintercept = 0, linetype = 2, size = 1) +
  geom_line(data = filter(SSB, retro != "retro_0"), 
            aes(x = Year, y = perc_diff, colour = retro, group = retro), size = 1) +
  scale_color_grey(guide = FALSE) +
  labs(x = NULL, y = "Percent difference\nfrom terminal year\n") +
  scale_y_continuous(limits = c(-50, 50)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels)-> retro2

cowplot::plot_grid(retro1, retro2, align = "hv", nrow = 2) -> retro_plot
retro_plot

ggsave(paste0(run_dir, "/retro_", MODEL_VERSION, ".png"), plot = retro_plot, dpi = 300, height = 6, width = 6, units = "in")

# Mohn's Rho ----

# Alaska Fisheries Science Center and Hurtado-Ferro et al. (2015) Mohn's rho
# https://www.afsc.noaa.gov/REFM/stocks/Plan_Team/2013/Sept/Retrospectives_2013_final3.pdf
# # mean over all peels (Estimate in peel year - reference estimate (current
# # year's estimate) / reference estimate)

bind_cols(SSB %>% 
            filter(Year == max(Year)) %>% 
            mutate(m = (spB - term_spB) / term_spB ) %>% 
            ungroup() %>% 
            summarize(mohns_rho = mean(m)),
          
          # Wood's Hole Rho (uses all years, not just the terminal year of the peel.
          # LeGault et al 2009 found no marked differences between the two rhos)
          SSB %>% 
            mutate(m = (spB - term_spB) / term_spB ) %>% 
            ungroup() %>% 
            summarize(woods_rho = mean(m)) ) -> rhos

write_csv(x = rhos, path = paste0(run_dir, "/retrospective_rhos.csv"))
# Interpretation of Mohn's Rho (Hurtado-Ferro et al 2015):
# Given that the variability of Mohn’s rdepends on life history, and that the
# statistic appears insensitive to F, we propose the following rule of thumb
# when determining whether a retrospective pattern should be addressed
# explicitly: values of Mohn’s r higher than 0.20 or lower than 20.15 for
# longer-lived species (upper and lower bounds of the 90% simulation intervals
# for the flatfish base case), or higher than 0.30 or lower than 20.22 for
# shorter-lived species (upper and lower bounds of the 90% simulation intervals
# for the sardine base case) should be cause for concern and taken as indicators
# of retrospective patterns. However, Mohn’s r values smaller than those
# proposed should not be taken as confirmation that a given assessment does not
# present a retrospective pattern, and the choice of 90% means that a “false
# positive” will arise 10% of the time.
