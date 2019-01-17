
# Sensitivity analysis on sigmaR (variability in recruitment deviations) and
# sigmaM (variability block-specific natural mortality deviations)

# Pseudocode 
# 1. Set up directories (create a subdirectory of HER called sensitivity_sigmas
# and copy the tpl to each run folder)
# 2. Select values of precision of recruitment deviations (precR = 1/(sigmaR^2))
# and sigma_M over which to run sensitivity analysis
# 3. Use a nested for loop that fixes values of precR over a range of sigmaM
# (and vice versa), while keeping track of convergence and model results

# Libraries/source files ----

library(R2admb) # run ADMB in R
library(tidyverse)
library(data.table)

source("R/create_ctl.R") # function that writes new ctl file
source("R/helper.R")

options(scipen = 999) # remove scientific notation

# Set up directories ----
main_dir <- getwd()#"S:\\Region1Shared-DCF\\Research\\Herring-Dive Fisheries\\Herring\\ADMB Rudd Workshop 2018\\AlaskaHerring"
proj_dir <- file.path(main_dir, "HER")
setwd(proj_dir)

files <- list.files(file.path(proj_dir)) # searchable file list

tpl <- file.path(proj_dir, files[which(grepl(".tpl", files))]) # find tpl
name <- strsplit(files[grepl("tpl",files)],".tpl")[[1]][1] # program name

dat <- file.path(proj_dir, files[which(grepl(".dat", files))]) # find dat(two)


# Compile tpl if it hasn't been compiled yet
exe <- paste0(name,".exe")
if (file.exists(exe) == FALSE) {
  setup_admb() # check admb path
  compile_admb(name, verbose = TRUE) # compile
} 

std <- paste0(name,".std") # for checking model convergence in loop

# Run analysis ----

# Variability of recruitment devs defined by precision = precR = 1/(sigmaR^2) or sigmaR
# = sqrt(1/precR). Smaller precisions result in large sigmaR and vice versa.
# Steve Martell fixed precR at 6.25 (sigmaR = 0.4)
# 

# Standard deviation of natural mortality devs. Steve Martell fixed sigmaM at
# 0.001
# 

# Values of sigmaR (precR) and sigmaM over which to run a sensitivity analysis:
sigmaR_vec <- seq(0.3, 1.7, by = 0.2)
precR_vec <- round(1/(sigmaR_vec^2), 3)
sigmaM_vec <- c(0.001, 0.005, 0.010, 0.050, 0.100, 0.500)

rundir <- file.path(proj_dir, "sensitivity_sigmas")
dir.create("sensitivity_sigmas", showWarnings = FALSE)

# check_converge <- list()
check_converge <- matrix(nrow = length(precR_vec), ncol = length(sigmaM_vec))

for(i in 1:length(precR_vec)){
  for(j in 1:length(sigmaM_vec)){

    sensdir <- file.path(rundir, paste0("precR_", precR_vec[i], "_sigmaM_", sigmaM_vec[j]))
    dir.create(sensdir, showWarnings = FALSE)
    setwd(sensdir)
    
    # Create new sitka_2.ctl file
    create_ctl(precR = precR_vec[i], sigmaM = sigmaM_vec[j])

    # copy TPL to run directory
    file.copy(from = tpl, to = sensdir, overwrite = TRUE)	
    
    # copy exe
    file.copy(from = file.path(proj_dir,exe), to = sensdir, overwrite = TRUE)	
    
    # copy data files
    file.copy(from = dat[1], to = sensdir, overwrite = TRUE)	
    file.copy(from = dat[2], to = sensdir, overwrite = TRUE)	
    
    # remove pre-existing std file if it exists
    if (file.exists(std)) {file.remove(std, showWarnings = FALSE)}
  
    # run admb, set seed at retro value
    run_admb(name, verbose = TRUE)
  
  # check for convergence (1 = converge, 0 = not converged)
  if (file.exists(std)) {
    check_converge[i,j] <- 1
  } else {check_converge[i,j] <- 0 }
  
  }
}


# Model selection ----

# Pseudocode for model selection:
# 1. Characterize and remove models that did not converge (both from check file
# where 0 = not converged, but also from max gradient).
# 2. Of remaining models, check model estimates of recruitment and natural mortality

check <- as.data.frame(check_converge, row.names = FALSE)
names(check) <- sigmaM_vec
check %>% 
  mutate(precR = precR_vec) %>% 
  gather("sigmaM", "converged?", -precR) %>% 
  mutate(fn = paste0("precR_", precR, "_sigmaM_", sigmaM)) -> check

# Save relevant files for later use 
setwd(rundir)

if(exists("check")) {
  write_csv(check, "sensitivity_analysis_results.csv")
  } else {
    check <- read_csv("sensitivity_analysis_results.csv") }

# 1. Characterize models that did not converge: Roughly 30% of the trials did
# not converge. These models were characterized by very high precision on
# recruitment (low sigmaR); none of the models with precR=11.1 (sigmaR=0.3)
# converged. Otherwise it was hard to detect a pattern with precR; a few from
# each precR value did not converge across the range of values explored.
# Interestingly model formulations with intermediate to high values of sigmaM
# were the least likely to converge (0.05, 0.1, 0.5), while models with low
# sigmaM (0.001) converged in all trials.
filter(check, `converged?` == 0) %>% arrange(precR)
filter(check, `converged?` == 0 & precR != max(precR)) %>% arrange(precR)
filter(check, `converged?` == 0 & precR != max(precR)) %>% arrange(sigmaM)

# Read results for converged models
filter(check, `converged?` == 1) -> check

diags_ls <- list() # store diagnostic output
derived_ls <- list() # store derived time series of interest
natmat_ls <- list() # time-varying natural mortality & survival
matsel_ls <- list() # time-varying maturity and selectivity

# NOTE- Year for output will some times indicate terminal year in a time block
for(i in 1:length(check$fn)){
  
  precR <- check$precR[i]
  sigmaM <- check$sigmaM[i]
  fn <- check$fn[i]
  
  sensdir <- file.path(rundir, fn)
  setwd(sensdir)
  
  D <- read_admb(name)
  
  # Diagnostics, key derived forecast quantities and parameter estimates
  diags_ls[[i]] <- data.frame(fn = fn,
                              precR = precR,
                              sigmaM = sigmaM,
                              nll = D[["fit"]]$nlogl,
                              maxgrad = D[["fit"]]$maxgrad,
                              nopar = D[["fit"]]$nopar,
                              ghl = D[["ghl"]],
                              fore_sb = D[["fore_sb"]] / 0.90718,
                              fore_matb = D[["fore_matb"]] / 0.90718,
                              rinit = exp(D[["theta"]][2]),
                              rbar = exp(D[["theta"]][3]))
  
  # Natural mortality and survival
  df <- data.frame(D[["year"]], 
                   D[["Mij"]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  
  natmat_ls[[i]] <- df %>% 
    select(Year, M = `3`) %>% # nat mat not age-specific, so just use one age
    mutate(survival = exp(-M),
           precR = precR,
           sigmaM = sigmaM,
           fn = fn) %>% 
    filter(Year %in% c(1998, 2014, 2017)) %>% # Predefined time blocks *FLAG* this is hard coded
    mutate(`Time blocks` = c("1980-1998", "1999-2014", "2015-2017"),
           `Time blocks` = factor(`Time blocks`, ordered = TRUE))
  
  # Maturity/selectivity
  df <- data.frame(D[["year"]], 
                   D[["mat"]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  df %>% mutate(param = "Maturity") -> df
  
  nyr <- D[["mod_nyr"]] - D[["mod_syr"]] + 1
  
  df2 <- data.frame(D[["year"]], 
                    D[["Sij"]][1:nyr, ]) # change number of rows
  colnames(df2) <- c("Year", paste(D[['sage']]:D[['nage']]))
  df2 %>% mutate(param = "Selectivity") %>% 
    bind_rows(df) -> df
  
  matsel_ls[[i]] <- df %>% 
    gather("Age", "proportion", -c(Year, param)) %>% 
    mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+")),
           # Normalize selectivity from 0 to 1
           proportion = ifelse(param == "Selectivity", (proportion - 0)/(max(proportion) - 0), 
                               proportion)) %>% 
    group_by(Age, param, proportion) %>% 
    mutate(min = min(Year),
           max = max(Year),
           `Time blocks` = paste0(min, "-", max),
           precR = precR,
           sigmaM = sigmaM,
           fn = fn) %>%
    ungroup() %>% 
    distinct(precR, sigmaM, fn, param, `Time blocks`, Year = max, Age, proportion)
  
  # Derived time series of biomass and ASA vs Ricker-predicted recruitment
    derived_ls[[i]] <- data.frame(Year = D[["year"]], 
                            spB = D[["sp_B"]] / 0.90718, # convert to short tons
                            matB = D[["mat_B"]] / 0.90718,
                            age3 = D[["Nij"]][1:38, 1] # ASA age-3 recruits
                            ) %>% 
      left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                      D[["mod_nyr"]] + 1, 1),
                           # SR is the current year age-3
                           SR = D[["recruits"]],
                           # ssb is the spawning biomass in year-3 that produced
                           # the current years recruits
                           ssb = D[["spawners"]] / 0.90718), by = "Year") %>% 
      mutate(precR = precR,
             sigmaM = sigmaM,
             fn = fn)
}

diags <- do.call(rbind, diags_ls)
derived <- do.call(rbind, derived_ls)
natmat <- do.call(rbind, natmat_ls)
matsel <- do.call(rbind, matsel_ls)

# Not converged! Remove these from all subsequent analyses
diags %>% filter(maxgrad > 0.001) %>% pull(fn) -> high_maxgrad
diags %>% filter(!fn %in% high_maxgrad) -> diags
derived %>% filter(!fn %in% high_maxgrad) -> derived
natmat %>% filter(!fn %in% high_maxgrad) -> natmat
matsel %>% filter(!fn %in% high_maxgrad) -> matsel

# Colour palettes ----

# For showing sigmaR gradient
palR <- colorRampPalette(c("#fdd0a2","#fa9fb5", "#4a1486"))
sigmaR_cols <- palR(length(sigmaR_vec))

# For showing sigmaM gradient
pal <- colorRampPalette(c("#edf8b1","#7fcdbb", "#2c7fb8"))
sigmaM_cols <- pal(length(sigmaM_vec))

# Diagnostics ----

setwd(proj_dir)

figdir <- file.path(proj_dir, "figs/sensitivity_sigmas")
dir.create("figs/sensitivity_sigmas", showWarnings = FALSE)

diags %>% mutate(sigmaR = round(sqrt(1/precR),1)) -> diags

diags$nopar # check that in fact all models have the same number of parameters

# Maximum gradient component
diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigmaR), size = maxgrad)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Maximum\ngradient",
             breaks = seq(min(diags$maxgrad), max(diags$maxgrad), length.out = 3),
             labels = round(seq(min(diags$maxgrad), max(diags$maxgrad), length.out = 3), 5)) +
  geom_point(data = filter(diags, maxgrad == min(maxgrad)),
             col = "green", size = 10, fill = NA, shape = 6) +
  geom_point(data = filter(diags, maxgrad == max(maxgrad)),
             col = "red", size = 10, fill = NA, shape = 2) +
  labs(x = "\nsigmaM", y = "sigmaR\n")

ggsave(paste0(figdir, "/maxgrad.png"), dpi = 300, height = 4, width = 5, units = "in")

# Objective function - nothing about this is surprisingly or even that
# informative. If estimated, these parameters would shrink to 0 because less
# variability is going to result in a better fit. The fact that this doesn't
# result in a huge spread of NLL values is a good sign.
diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigmaR), size = nll)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(trans = "reverse",
             range = c(0, 4),
             name = "Negative\nlog\nlikelihood",
             breaks = seq(min(diags$nll), max(diags$nll), length.out = 3),
             labels = round(seq(min(diags$nll), max(diags$nll), length.out = 3), 0)) +
  labs(x = "\nsigmaM", y = "sigmaR\n")

ggsave(paste0(figdir, "/nll.png"), dpi = 300, height = 4, width = 5, units = "in")

# Initial mean numbers-at-age (rinits) - as might be expected, the higher sigmaM is,
# the larger initial recruitment is. Increasing sigmaR also increases rinit, but
# sigmaM is far more influential.
diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigmaR), size = rinit)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Initial mean\nnumbers-at-age\n(millions)",
             breaks = seq(min(diags$rinit), max(diags$rinit), length.out = 5),
             labels = round(seq(min(diags$rinit), max(diags$rinit), length.out = 5), 0)) +
  labs(x = "\nsigmaM", y = "sigmaR\n")

ggsave(paste0(figdir, "/rinit.png"), dpi = 300, height = 4, width = 5, units = "in")

# Age-3 mean recruitment (rbar) - Exact same results and rationale as initial
# mean numbers-at-age
diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigmaR), size = rbar)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Age-3\nmean\nrecruitment\n(millions)",
             breaks = seq(min(diags$rbar), max(diags$rbar), length.out = 5),
             labels = round(seq(min(diags$rbar), max(diags$rbar), length.out = 5), 0)) +
  labs(x = "\nsigmaM", y = "sigmaR\n")

ggsave(paste0(figdir, "/rbar.png"), dpi = 300, height = 4, width = 5, units = "in")

# GHL ----

# Large sigmaM's result in drastically smaller GHLs overall. The minimum,
# maximum, and approximate median values are shown in red triangle, green
# upsidedown triangle, and orange square.

diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigmaR), size = ghl)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "GHL (t)",
             breaks = round(seq(min(diags$ghl), max(diags$ghl), length.out = 5), 0),
             labels = comma) + 
  geom_point(data = filter(diags, ghl == max(ghl)),
             col = "red", size = 10, fill = NA, shape = 2) +
  geom_point(data = filter(diags, ghl == min(ghl)),
             col = "green", size = 10, fill = NA, shape = 6) +
  geom_point(data = slice(diags, which.min(abs(ghl-median(ghl)))),
             col = "orange", size = 10, fill = NA, shape = 0) +
  labs(x = "\nsigmaM", y = "sigmaR\n")

ggsave(paste0(figdir, "/ghl.png"), dpi = 300, height = 4, width = 5, units = "in")

# Natural mortality ----

natmat %>% mutate(sigmaR = round(sqrt(1/precR), 1)) -> natmat

# increasing sigmaM has dramatic effects on time-varying M
ggplot(natmat, aes(x = `Time blocks`, y = survival, 
                   group = sigmaM, colour = sigmaM)) +
  geom_point(size = 1.3) + 
  geom_line(size = 1.3) + 
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~sigmaR) +
  labs(colour = "sigmaM", y = "Survival")

ggsave(paste0(figdir, "/survival_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

# sigmaR has an affect on M at values of sigmaM < 0.01 (low sigmaR results in
# higher survival), but the affect is swamped at high values of sigmaM
ggplot(natmat, aes(x = `Time blocks`, y = survival, 
                   group = factor(sigmaR), colour = factor(sigmaR))) +
  geom_line(size = 1.3) + 
  geom_point(size = 1.3) + 
  scale_colour_manual(values = sigmaR_cols) +
  facet_wrap(~sigmaM) +
  labs(colour = "sigmaR", y = "Survival")

ggsave(paste0(figdir, "/survival_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# Maturity and selectivity ----
matsel %>% 
  mutate(sigmaR = round(sqrt(1/precR), 1)) %>% 
  group_by(param) %>% 
  mutate(combos = paste0("sigmaM: ", sigmaM, " (", `Time blocks`, ")"),
         combos_R = paste0("sigmaR: ", sigmaR, " (", `Time blocks`, ")")) -> matsel

# Trends appear constant across sigmaR. Increasing sigmaM appears to result in
# the maturity curves shifting to the right (maturing later)
matsel %>% filter(param == "Maturity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = sigmaM, 
                group = combos), size = 1) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~sigmaR) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Maturity\n", colour = "sigmaM")

ggsave(paste0(figdir, "/maturity_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

# Confirms that there is very slight variation in maturity across sigmaR for
# small values of sigmaM), where maturation occurs slightly later with increased
# sigmaR in the first time block (1980-2014).
matsel %>% filter(param == "Maturity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = factor(sigmaR), 
                group = combos_R), size = 1) +
  scale_colour_manual(values = sigmaR_cols) +
  facet_wrap(~sigmaM) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Maturity\n", colour = "sigmaR")

ggsave(paste0(figdir, "/maturity_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# sigmaM and sigmaR do not influence selectivity estimates
matsel %>% filter(param == "Selectivity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = sigmaM, 
                group = combos), size = 1) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~sigmaR) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Selecitivity\n", colour = "sigmaM")

ggsave(paste0(figdir, "/selectivity_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

matsel %>% filter(param == "Selectivity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = factor(sigmaR), 
                group = combos_R), size = 1) +
  scale_colour_manual(values = sigmaR_cols) +
  facet_wrap(~sigmaM) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Selecitivity\n", colour = "sigmaR")

ggsave(paste0(figdir, "/selectivity_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# Derived time series ----
derived %>% mutate(sigmaR = round(sqrt(1/precR), 1)) -> derived

# Separate forecast quantities into separate df
diags %>% 
  select(fn, sigmaR, sigmaM, fore_sb, fore_matb) %>% 
  mutate(Year = max(derived$Year) + 1) -> forec

# SSB facetted by sigmaM - shows variability in ssb is attributed to sigmaM not
# sigmaR

# Facet by sigmaR (colours are sigmaM)
M <- derived %>% 
  ggplot(aes(x = Year, colour = sigmaM, group = sigmaM)) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~sigmaR) +
  scale_y_continuous(label = scales::comma)

# Facet by sigmaM (colours are sigmaR)
R <- derived %>% 
  ggplot(aes(x = Year, colour = factor(sigmaR), group = factor(sigmaR))) +
  scale_colour_manual(values = sigmaR_cols) +
  facet_wrap(~sigmaM) +
  scale_y_continuous(label = scales::comma)

# Spawning/mature biomass -----

# sigmaM strongly influences estimates; low sigmaM results in lower estimates at
# beginning of ts, then much higher during the 2000s.
M + geom_line(aes(y = spB), size = 1) +
  labs(x = "\nYear", y = "Spawning biomass (t)\n", colour = "sigmaM") +
  geom_point(data = forec, aes(x = Year, y = fore_sb, colour = sigmaM, group = sigmaM),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/spbiomass_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

# sigmaR has no influence on estimates
R + geom_line(aes(y = spB), size = 1) +
  labs(x = "\nYear", y = "Spawning biomass (t)\n", colour = "sigmaR")

ggsave(paste0(figdir, "/spbiomass_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# Mature biomass: same results as spawning biomass
M + geom_line(aes(y = matB), size = 1) +
  labs(x = "\nYear", y = "Mature biomass (t)\n", colour = "sigmaM")

ggsave(paste0(figdir, "/matbiomass_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

R + geom_line(aes(y = matB), size = 1) +
  labs(x = "\nYear", y = "Mature biomass (t)\n", colour = "sigmaR")

ggsave(paste0(figdir, "/matbiomass_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# Recruitment ----

# ASA-estimated age-3 abundance: 

# sigmaM appears to have a slight affect on estimates but trends are the same.
M + geom_line(aes(y = age3), size = 1) +
  labs(x = "\nYear", y = "ASA age-3 abundance (millions)\n", colour = "sigmaM")

# sigmaR has no affect on estimates
R + geom_line(aes(y = age3), size = 1) +
  labs(x = "\nYear", y = "ASA age-3 abundance (millions)\n", colour = "sigmaR")

# Ricker-estimated age-3 abundance:

# Very different trends in first part of ts, similar after 2000-ish.
M + geom_line(aes(y = SR), size = 1) +
  labs(x = "\nYear", y = "Ricker age-3 abundance (millions)\n", colour = "sigmaM")

# sigmaR highly influential. As expected, higher sigmaR results in more
# variability over the ts. Trends over time are the same for a given sigmaM
R + geom_line(aes(y = SR), size = 1) +
  labs(x = "\nYear", y = "Ricker age-3 abundance (millions)\n", colour = "sigmaR")

head(diags)

