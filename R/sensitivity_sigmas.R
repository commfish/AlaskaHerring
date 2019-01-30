
# Sensitivity analysis on precR/sigmaR (precision/variability in stock-recruitment
# relationship), sigma_rdevs (variability in recruitment deviations), and
# sigmaM (variability block-specific natural mortality deviations)

# Definitions


# precR - precision of stock-recruitment relationship (SM value:
# precR=6.25, sigmaR=0.4, JS changed to: precR=4.0m sigmaR=0.5 b/c more trials converge across range of
# sigmaM)

# sigma_rdevs (and f_sigma_rdevs) - variability of age-3 recruitment deviations
# (rbar_devs) and initial numbers-at-age (rinit_devs) during all phases of
# estimation vs. the final phase (SM value: sigma_rdevs=1.0, f_sigma_rdevs=5.0)

# sigmaM - variability of natural mortality deviations

# Pseudocode 
# 1. Set up directories (create a subdirectory of HER called sensitivity_sigmas
# and copy the tpl to each run folder)
# 2. Select values of precision of recruitment deviations (precR = 1/(sigmaR^2))
# and sigma_M over which to run sensitivity analysis
# 3. Use a nested for loop that fixes values of precR over a range of sigmaM
# (and vice versa), while keeping track of convergence and model results
# 4. Select a sigmaR (chose 0.5 over 0.4 because it improved convergence)
# 5. Select values of sigma_rdevs/f_sigma_rdevs such that we capture low,
# intermediate and high values of sigma_rdevs crossed with fixed versus phased
# f_sigma_rdevs crossed with low and high variability of f_sigma_rdevs
# 6. Use a nested for loop over to examine relationship between sigmaM and the
# sigma_rdevs

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

# Standard deviation of natural mortality devs. Steve Martell fixed sigmaM at
# 0.001

# Values of sigmaR (precR) and sigmaM over which to run a sensitivity analysis:
sigmaR_vec <- seq(0.3, 1.6, by = 0.2)
sigmaR_vec <- c(0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.4, 2.0)
sigmaR_vec <- c(0.4, 0.5, 0.6)
precR_vec <- round(1/(sigmaR_vec^2), 3)
sigmaM_vec <- c(0.001, 0.010, 0.025, 0.050, 0.075, 0.100, 0.500, 1.00)

# For looking at sigmaR and sigmaM
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
                              ghl = D[["ghl"]], # Already in short tons
                              fore_sb = D[["fore_sb"]] / 0.90718, # convert to short tons
                              fore_matb = D[["fore_matb"]] / 0.90718, # convert to short tons
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
palR <- colorRampPalette(c("#bfd3e6","#8c6bb1", "#4d004b"))
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
  geom_line(size = 1) + 
  geom_point(size = 1.3) + 
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~sigmaR) +
  labs(colour = "sigmaM", y = "Survival")

ggsave(paste0(figdir, "/survival_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

# sigmaR has an effect on M at values of sigmaM < 0.01 (low sigmaR results in
# higher survival), but the effect is swamped at high values of sigmaM
ggplot(natmat, aes(x = `Time blocks`, y = survival, 
                   group = factor(sigmaR), colour = factor(sigmaR))) +
  geom_line(size = 1) + 
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

ggsave(paste0(figdir, "/maturity_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

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

ggsave(paste0(figdir, "/selectivity_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

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
derived %>%
  mutate(sigmaR = round(sqrt(1/precR), 1)) -> derived

# Separate forecast quantities into separate df - cannot use the model-supplied
# fore_sb b/c the assumptions are not correct
diags %>% 
  select(fn, sigmaR, sigmaM, fore_matb, ghl) %>% 
  mutate(Year = max(derived$Year) + 1,
         fore_sb = fore_matb - ghl) -> forec

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
M + geom_line(aes(y = spB), size = 0.8) +
  labs(x = "\nYear", y = "Spawning biomass (t)\n", colour = "sigmaM") +
  geom_point(data = forec, aes(x = Year, y = fore_sb, colour = sigmaM, group = sigmaM),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/spbiomass_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

# sigmaR has no influence on estimates
R + geom_line(aes(y = spB), size = 0.8) +
  labs(x = "\nYear", y = "Spawning biomass (t)\n", colour = "sigmaR") +
  geom_point(data = forec, aes(x = Year, y = fore_sb, colour = factor(sigmaR), group = factor(sigmaR)),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/spbiomass_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# Mature biomass: same results as spawning biomass
M + geom_line(aes(y = matB), size = 0.8) +
  labs(x = "\nYear", y = "Mature biomass (t)\n", colour = "sigmaM") +
  geom_point(data = forec, aes(x = Year, y = fore_matb, colour = sigmaM, group = sigmaM),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/matbiomass_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

R + geom_line(aes(y = matB), size = 0.8) +
  labs(x = "\nYear", y = "Mature biomass (t)\n", colour = "sigmaR") +
  geom_point(data = forec, aes(x = Year, y = fore_matb, colour = factor(sigmaR), group = factor(sigmaR)),
             shape = "*", size = 4) 

ggsave(paste0(figdir, "/matbiomass_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# Recruitment ----

# ASA-estimated age-3 abundance: 

# sigmaM appears to have a slight affect on estimates but trends are the same.
M + geom_line(aes(y = age3), size = 0.8) +
  labs(x = "\nYear", y = "ASA age-3 abundance (millions)\n", colour = "sigmaM")

ggsave(paste0(figdir, "/ASA_age3_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

# sigmaR has no affect on estimates
R + geom_line(aes(y = age3), size = 0.8) +
  labs(x = "\nYear", y = "ASA age-3 abundance (millions)\n", colour = "sigmaR")

ggsave(paste0(figdir, "/ASA_age3_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# Ricker-estimated age-3 abundance:

# Very different trends in first part of ts, similar after 2000-ish.
M + geom_line(aes(y = SR), size = 0.8) +
  labs(x = "\nYear", y = "Ricker age-3 abundance (millions)\n", colour = "sigmaM")

ggsave(paste0(figdir, "/Ricker_age3_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

# sigmaR highly influential. As expected, higher sigmaR results in more
# variability over the ts. Trends over time are the same for a given sigmaM
R + geom_line(aes(y = SR), size = 0.8) +
  labs(x = "\nYear", y = "Ricker age-3 abundance (millions)\n", colour = "sigmaR")

ggsave(paste0(figdir, "/Ricker_age3_sigmaR.png"), dpi = 300, height = 4, width = 8, units = "in")

# Recommondations decrease precR to 4.0 from 6.25 b/c it increases convergence
# probability across sigmaM increase sigmaM to a range ~0.075 b/c it appears to
# provide an intermediate level of variability for M while still producing
# reasonable estimates (especially in the last time period that is very short)

# Analysis II ----

# Sensitivity analysis of sigma_rdevs/f_sigma_rdevs and sigmaM

# Scheme reviewed and approved by SD 2019-01-25. Designed to capture low,
# medium, high variability crossed with different phasing crossed with low and
# high variability for the final phase. Book-ended by what I deemed biologically
# reasonable for herring.
sigma_rdevs_df <- data.frame(description = c("Low initial variability in rec devs, phasing fixed",
                                             "Low initial variability in rec devs, final phase low variability",
                                             "Low initial variability in rec devs, final phase high variability",
                                             "Intermediate initial variability in rec devs, phasing fixed",
                                             "Intermediate initial variability in rec devs, final phase low",
                                             "Intermediate initial variability in rec devs, final phase high",
                                             "High initial variability in rec devs, phasing fixed",
                                             "High initial variability in rec devs, final phase low",
                                             "High initial variability in rec devs, final phase high"),
                             sigma_rdevs = c(rep(0.8, 3),
                                             rep(1.0, 3),
                                             rep(1.6, 3)),
                             f_sigma_rdevs = c(0.8, 2, 5,
                                               1.0, 2, 5,
                                               1.6, 2, 5))
sigma_rdevs_df %>% 
  mutate(combos = paste(sigma_rdevs, f_sigma_rdevs, sep = "_")) -> sigma_rdevs_df

sigmaM_vec <- c(0.001, 0.010, 0.020, 0.030, 0.040,
                0.050, 0.060, 0.070, 0.080, 0.090, 
                0.100, 0.200, 0.300, 0.400, 0.500)

# For looking at sigmaR and sigmaM
rundir <- file.path(proj_dir, "sensitivity_rdevs_sigmaM")
dir.create("sensitivity_rdevs_sigmaM", showWarnings = FALSE)

check_converge <- matrix(nrow = nrow(sigma_rdevs_df), ncol = length(sigmaM_vec))

for(i in 1:nrow(sigma_rdevs_df)){
  for(j in 1:length(sigmaM_vec)){
    
    sensdir <- file.path(rundir, paste0("rdevs_", sigma_rdevs_df$combos[i], 
                                        "_sigmaM_", sigmaM_vec[j]))
    dir.create(sensdir, showWarnings = FALSE)
    setwd(sensdir)
    
    # Create new sitka_2.ctl file
    create_ctl(precR = 4.0, sigmaM = sigmaM_vec[j],
               sigma_rdevs = sigma_rdevs_df$sigma_rdevs[i],
               f_sigma_rdevs = sigma_rdevs_df$f_sigma_rdevs[i])
    
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


check <- as.data.frame(check_converge, row.names = FALSE)
names(check) <- sigmaM_vec
check %>% 
  mutate(rdevs_combos = sigma_rdevs_df$combos,
         sigma_rdevs = sigma_rdevs_df$sigma_rdevs,
         f_sigma_rdevs = sigma_rdevs_df$f_sigma_rdevs) %>% 
  gather("sigmaM", "converged?", -c(sigma_rdevs, f_sigma_rdevs, rdevs_combos)) %>% 
  mutate(fn = paste0("rdevs_", rdevs_combos, 
                     "_sigmaM_", sigmaM)) -> check

# II. Model selection -----

#Characteristics of convergence:

# All trials with sigmaM <= 0.03 converged.
# sigmaM = 0.04 no trials converged.
# sigmaM = 0.05 all trials converged.
# sigmaM >= 0.06 mixed convergence patterns. Best convergence performance at
# sigmaM=c(0.08,0.4,0.5) with 2/3 trials converged

# Convergence issues in midrange of sigmaM (which is what we think might be most
# biologically reasonable). Likely due to constrained parameterization of M and
# poorly estimated deviation in the final time block.

# No pattern was revealed for convergence by sigma_rdev or f_sigma_rdev. All
# combinations result in roughly 2/3 convergence across sigmaM values.

table(check$f_sigma_rdevs, check$`converged?`)
table(check$sigma_rdevs, check$`converged?`)
table(check$sigmaM, check$`converged?`)
table(check$rdevs_combos, check$`converged?`)

# Save relevant files for later use 
setwd(rundir)

if(exists("check")) {
  write_csv(check, "sensitivity_rdevs_sigmaM_results.csv")
} else {
  check <- read_csv("sensitivity_rdevs_sigmaM_results.csv") }

# Read results for converged models
filter(check, `converged?` == 1) -> check

diags_ls <- list() # store diagnostic output
derived_ls <- list() # store derived time series of interest
natmat_ls <- list() # time-varying natural mortality & survival
matsel_ls <- list() # time-varying maturity and selectivity

# NOTE- Year for output will some times indicate terminal year in a time block
for(i in 1:length(check$fn)){
  
  sigma_rdevs <- check$sigma_rdevs[i]
  f_sigma_rdevs <- check$f_sigma_rdevs[i]
  rdevs_combo <- check$rdevs_combos[i]
  sigmaM <- check$sigmaM[i]
  fn <- check$fn[i]
  
  sensdir <- file.path(rundir, fn)
  setwd(sensdir)
  
  D <- read_admb(name)
  
  # Diagnostics, key derived forecast quantities and parameter estimates
  diags_ls[[i]] <- data.frame(fn = fn,
                              sigma_rdevs = sigma_rdevs,
                              f_sigma_rdevs = f_sigma_rdevs,
                              rdevs_combo = rdevs_combo,
                              sigmaM = sigmaM,
                              nll = D[["fit"]]$nlogl,
                              maxgrad = D[["fit"]]$maxgrad,
                              nopar = D[["fit"]]$nopar,
                              ghl = D[["ghl"]], # Already in short tons
                              fore_sb = D[["fore_sb"]] / 0.90718, # convert to short tons
                              fore_matb = D[["fore_matb"]] / 0.90718, # convert to short tons
                              rinit = exp(D[["theta"]][2]),
                              rbar = exp(D[["theta"]][3]))
  
  # Natural mortality and survival
  df <- data.frame(D[["year"]], 
                   D[["Mij"]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  
  natmat_ls[[i]] <- df %>% 
    select(Year, M = `3`) %>% # nat mat not age-specific, so just use one age
    mutate(survival = exp(-M),
           sigma_rdevs = sigma_rdevs,
           f_sigma_rdevs = f_sigma_rdevs,
           rdevs_combo = rdevs_combo,
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
           sigma_rdevs = sigma_rdevs,
           f_sigma_rdevs = f_sigma_rdevs,
           rdevs_combo = rdevs_combo,
           sigmaM = sigmaM,
           fn = fn) %>%
    ungroup() %>% 
    distinct(sigma_rdevs, f_sigma_rdevs, rdevs_combo, sigmaM, fn, param, 
             `Time blocks`, Year = max, Age, proportion)
  
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
    mutate(sigma_rdevs = sigma_rdevs,
           f_sigma_rdevs = f_sigma_rdevs,
           rdevs_combo = rdevs_combo,
           sigmaM = sigmaM,
           fn = fn)
}

diags <- do.call(rbind, diags_ls)
derived <- do.call(rbind, derived_ls)
natmat <- do.call(rbind, natmat_ls)
matsel <- do.call(rbind, matsel_ls)

# Not converged! Remove these from all subsequent analyses. There were 4 of
# these (3 sigmaM=0.1 and 1 sigmaM=0.3)
diags %>% filter(maxgrad > 0.001) %>% pull(fn) -> high_maxgrad

# Clean up function to get rid of any mods with high max gradient components and
# create new factor based on phasing treatment of sigma_rdevs for sensitivity
# analysis results
cleanup <- function(df) {
  df %>% 
    filter(!fn %in% high_maxgrad) %>% 
    mutate(`Phasing of sigma_rdevs` = derivedFactor(
      "Constant sigma_rdevs" = sigma_rdevs == f_sigma_rdevs,
      "Final sigma_rdevs = 2.0" = f_sigma_rdevs == 2.0,
      "Final sigma_rdevs = 5.0" = f_sigma_rdevs == 5.0,
      .default = NA)) %>% 
    mutate(rdevs_combo = factor(rdevs_combo,
                                levels = c("0.8_0.8", "0.8_2", "0.8_5",
                                           "1_1", "1_2", "1_5",
                                           "1.6_1.6", "1.6_2", "1.6_5"),
                                labels = c("0.8_0.8", "0.8_2.0", "0.8_5.0",
                                           "1.0_1.0", "1.0_2.0", "1.0_5.0",
                                           "1.6_1.6", "1.6_2.0", "1.6_5.0"),
                                ordered = TRUE),
           sigma_rdevs = factor(sigma_rdevs,
                                levels = c("0.8", "1", "1.6"),
                                labels = c("0.8", "1.0", "1.6"),
                                ordered = TRUE))
}

cleanup(diags) -> diags
cleanup(derived) -> derived
cleanup(natmat) -> natmat
cleanup(matsel) -> matsel

# II. Diagnostics ----

setwd(proj_dir)

figdir <- file.path(proj_dir, "figs/sensitivity_rdevs_sigmaM")
dir.create("figs/sensitivity_rdevs_sigmaM", showWarnings = FALSE)

diags$nopar # check that in fact all models have the same number of parameters

# Maximum gradient component - the models that increase sigma_rdevs in the final
# phase have the best max gradient component
diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigma_rdevs), size = maxgrad)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Maximum\ngradient",
             breaks = seq(min(diags$maxgrad), max(diags$maxgrad), length.out = 3),
             labels = round(seq(min(diags$maxgrad), max(diags$maxgrad), length.out = 3), 5)) +
  geom_point(data = filter(diags, maxgrad == min(maxgrad)),
             col = "green", size = 10, fill = NA, shape = 6) +
  geom_point(data = filter(diags, maxgrad == max(maxgrad)),
             col = "red", size = 10, fill = NA, shape = 2) +
  facet_wrap(~`Phasing of sigma_rdevs`, ncol = 1) + 
  labs(x = "\nsigmaM", y = "sigma_rdevs\n")

ggsave(paste0(figdir, "/maxgrad.png"), dpi = 300, height = 7, width = 6, units = "in")

# Objective function - models with low sigmaM and models that allow an increase
# of sigma_rdev to 5.0 in the final phase have the smallest obj fun
diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigma_rdevs), size = nll)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(trans = "reverse",
             range = c(0, 4),
             name = "Negative\nlog\nlikelihood",
             breaks = seq(min(diags$nll), max(diags$nll), length.out = 3),
             labels = round(seq(min(diags$nll), max(diags$nll), length.out = 3), 0)) +
  facet_wrap(~`Phasing of sigma_rdevs`, ncol = 1) + 
  labs(x = "\nsigmaM", y = "sigma_rdevs\n")

ggsave(paste0(figdir, "/nll.png"), dpi = 300, height = 7, width = 6, units = "in")

# Initial mean numbers-at-age (rinits) - as might be expected, the higher sigmaM
# is, the larger initial recruitment is. Increasing sigma_rdevs does not appear
# to influence rinit. Changing the phasing for sigma_rdevs has some effect, but
# it appears small (keeping sigma_rdevs constant across phases results in larger
# rinit). sigmaM is far more influential.
diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigma_rdevs), size = rinit)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Initial mean\nnumbers-at-age\n(millions)",
             breaks = seq(min(diags$rinit), max(diags$rinit), length.out = 5),
             labels = round(seq(min(diags$rinit), max(diags$rinit), length.out = 5), 0)) +
  geom_point(data = filter(diags, rinit == max(rinit)),
             col = "red", size = 10, fill = NA, shape = 2) +
  geom_point(data = filter(diags, rinit == min(rinit)),
             col = "green", size = 10, fill = NA, shape = 6) +
  geom_point(data = slice(diags, which.min(abs(rinit-median(rinit)))),
             col = "orange", size = 10, fill = NA, shape = 0) +
  facet_wrap(~`Phasing of sigma_rdevs`, ncol = 1) + 
  labs(x = "\nsigmaM", y = "sigma_rdevs\n")

ggsave(paste0(figdir, "/rinit.png"), dpi = 300, height = 7, width = 6, units = "in")

# Age-3 mean recruitment (rbar) - Exact same results and rationale as initial
# mean numbers-at-age
diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigma_rdevs), size = rbar)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Age-3 mean\nrecruitment\n(millions)",
             breaks = seq(min(diags$rbar), max(diags$rbar), length.out = 5),
             labels = round(seq(min(diags$rbar), max(diags$rbar), length.out = 5), 0)) +
  geom_point(data = filter(diags, rbar == max(rbar)),
             col = "red", size = 10, fill = NA, shape = 2) +
  geom_point(data = filter(diags, rbar == min(rbar)),
             col = "green", size = 10, fill = NA, shape = 6) +
  geom_point(data = slice(diags, which.min(abs(rbar-median(rbar)))),
             col = "orange", size = 10, fill = NA, shape = 0) +
  facet_wrap(~`Phasing of sigma_rdevs`, ncol = 1) + 
  labs(x = "\nsigmaM", y = "sigma_rdevs\n")

ggsave(paste0(figdir, "/rbar.png"), dpi = 300, height = 7, width = 6, units = "in")

# II. GHL ----

# Large sigmaM's result in smaller GHLs overall. The largest GHLs occur with the
# smallest sigmaM and by increasing sigma_rdevs to 5.0 in the final phase (this
# maximized GHL regardless of the initial value for sigma_rdevs). The minimum,
# maximum, and approximate median values are shown in red triangle, green
# upsidedown triangle, and orange square.

diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigma_rdevs), size = ghl)) +
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
  facet_wrap(~`Phasing of sigma_rdevs`, ncol = 1) + 
  labs(x = "\nsigmaM", y = "sigma_rdevs\n")

ggsave(paste0(figdir, "/ghl.png"), dpi = 300, height = 7, width = 6, units = "in")

diags %>% 
  ggplot(aes(x = sigmaM, y = factor(sigma_rdevs))) +
  geom_tile(aes(fill = ghl)) +
  geom_text(aes(label = prettyNum(round(ghl, 0), big.mark=",")), 
            size = 2.5) +
  scale_fill_gradient2(
             name = "GHL (t)",
             labels = comma,
             low = muted("forestgreen"), 
             mid = "yellow", 
             high = "darkred", 
             midpoint = mean(diags$ghl)) + 
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        plot.title = element_text(size=20,face="bold")) +
  facet_wrap(~`Phasing of sigma_rdevs`, ncol = 1) + 
  labs(x = "\nsigmaM", y = "sigma_rdevs\n")

ggsave(paste0(figdir, "/ghl2.png"), dpi = 300, height = 7, width = 7, units = "in")

# II. Colour palettes ----

# Remove of the sigmaM values so the results are a little easier to digest sigmaM=0.2 and 0.4 for remaining visualizations so it's easier to read
derived %>% filter(!sigmaM %in% c(0.06,0.3,0.4,0.5)) -> derived
natmat %>% filter(!sigmaM %in% c(0.06,0.3,0.4,0.5)) -> natmat
matsel %>% filter(!sigmaM %in% c(0.06,0.3,0.4,0.5)) -> matsel
diags %>% filter(!sigmaM %in% c(0.06,0.3,0.4,0.5)) -> diags

# For showing sigmaR gradient
palR <- colorRampPalette(c("#ffeda0","#feb24c", "#f03b20"))
rdev_cols <- palR(length(unique(diags$rdevs_combo)))

# For showing sigmaM gradient
pal <- colorRampPalette(c("#edf8b1","#7fcdbb", "#2c7fb8"))
sigmaM_cols <- pal(length(unique(derived$sigmaM)))

# II. Natural mortality ----

# Increasing sigmaM has dramatic effects on time-varying M. It appears that the
# deviation in the final time block may be poorly estimated. These results are
# very similar to what was seen in the sensitivity analysis with sigmaR
ggplot(natmat, aes(x = `Time blocks`, y = survival, 
                   group = sigmaM, colour = sigmaM)) +
  geom_line(size = 1) + 
  geom_point(size = 1.3) + 
  scale_colour_manual(values = sigmaM_cols) +
  facet_grid(`Phasing of sigma_rdevs` ~ sigma_rdevs) +
  labs(colour = "sigmaM", y = "Survival")

ggsave(paste0(figdir, "/survival_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

# It appears sigma_rdevs has a very small effect on M at values of sigmaM < 0.01
# (high sigma_rdevs results in a slightly higher survival), but the effect goes
# away with increasing sigmaM

ggplot(natmat, aes(x = `Time blocks`, y = survival, 
                   group = factor(rdevs_combo), 
                   colour = factor(rdevs_combo))) +
  geom_line(size = 1) + 
  geom_point(size = 1.3) + 
  scale_colour_manual(values = rdev_cols) +
  facet_wrap(~sigmaM) +
  labs(colour = "sigma_rdevs", y = "Survival")

ggsave(paste0(figdir, "/survival_sigma_rdevs.png"), dpi = 300, height = 6, width = 8, units = "in")

# II: Maturity and selectivity ----
matsel %>% 
  group_by(param) %>% 
  mutate(combos = paste0("sigmaM: ", sigmaM, " (", `Time blocks`, ")"),
         combos_rdevs = paste0("sigma_rdevs: ", rdevs_combo, " (", `Time blocks`, ")")) -> matsel

# Trends appear constant across sigma_rdevs. Increasing sigmaM appears to result
# in the maturity curves shifting to the right (maturing later).
matsel %>% filter(param == "Maturity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = sigmaM, 
                group = combos), size = 1) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_grid(`Phasing of sigma_rdevs` ~ sigma_rdevs) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 3.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Maturity\n", colour = "sigmaM")

ggsave(paste0(figdir, "/maturity_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

# sigma_rdevs do not influence maturity estimates
matsel %>% filter(param == "Maturity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = rdevs_combo, 
                group = combos_rdevs), size = 1) +
  scale_colour_manual(values = rdev_cols) +
  facet_wrap(~sigmaM) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 3.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Maturity\n", colour = "sigma_rdevs")

ggsave(paste0(figdir, "/maturity_sigma_rdevs.png"), dpi = 300, height = 6, width = 8, units = "in")

# sigmaM and sigma_rdevs do not influence selectivity estimates
matsel %>% filter(param == "Selectivity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = sigmaM, 
                group = combos), size = 1) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_grid(`Phasing of sigma_rdevs` ~ sigma_rdevs) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 3.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Selecitivity\n", colour = "sigmaM")

ggsave(paste0(figdir, "/selectivity_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

matsel %>% filter(param == "Selectivity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = rdevs_combo, 
                group = combos_rdevs), size = 1) +
  scale_colour_manual(values = rdev_cols) +
  facet_wrap(~sigmaM) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 3.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Selecitivity\n", colour = "sigmaR")

ggsave(paste0(figdir, "/selectivity_sigma_rdevs.png"), dpi = 300, height = 6, width = 8, units = "in")

# II. Derived time series ----

# Facet by sigma_rdevs (colours are sigmaM)
M <- derived %>% 
  ggplot(aes(x = Year, colour = sigmaM, group = sigmaM)) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_grid(`Phasing of sigma_rdevs` ~ sigma_rdevs) +
  scale_y_continuous(label = scales::comma)

# Facet by sigmaM (colours are sigmaR)
R <- derived %>% 
  ggplot(aes(x = Year, colour = rdevs_combo, group = rdevs_combo)) +
  scale_colour_manual(values = rdev_cols) +
  facet_wrap(~sigmaM) +
  scale_y_continuous(label = scales::comma)

# Forecast values
diags %>% 
  select(fn, sigma_rdevs, rdevs_combo, `Phasing of sigma_rdevs`,
         sigmaM, fore_matb, ghl) %>% 
  mutate(Year = max(derived$Year) + 1,
         fore_sb = fore_matb - ghl) -> forec

# II. Spawning/mature biomass -----

# sigmaM strongly influences estimates; low sigmaM results in lower estimates at
# beginning of the time series, then much higher during the 2000s. Overall,
# higher values of sigmaM result in less variability in the time series.
# Noteworthy is that survival is estimated to be very low in the final time
# block with causes the forecast to descrease with increasing sigmaM
M + geom_line(aes(y = spB), size = 0.8) +
  labs(x = "\nYear", y = "Spawning biomass (t)\n", colour = "sigmaM") +
  geom_point(data = forec, aes(x = Year, y = fore_sb, colour = sigmaM, group = sigmaM),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/spbiomass_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

# sigma_rdevs has no influence on estimates
R + geom_line(aes(y = spB), size = 0.8) +
  labs(x = "\nYear", y = "Spawning biomass (t)\n", colour = "sigma_rdevs") +
  geom_point(data = forec, aes(x = Year, y = fore_sb, colour = rdevs_combo, group = rdevs_combo),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/spbiomass_sigma_rdevs.png"), dpi = 300, height = 6, width = 8, units = "in")

# Mature biomass: same results as spawning biomass
M + geom_line(aes(y = matB), size = 0.8) +
  labs(x = "\nYear", y = "Mature biomass (t)\n", colour = "sigmaM") +
  geom_point(data = forec, aes(x = Year, y = fore_matb, colour = sigmaM, group = sigmaM),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/matbiomass_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

R + geom_line(aes(y = matB), size = 0.8) +
  labs(x = "\nYear", y = "Mature biomass (t)\n", colour = "sigma_rdevs") +
  geom_point(data = forec, aes(x = Year, y = fore_matb, colour = rdevs_combo, group = rdevs_combo),
             shape = "*", size = 4) 

ggsave(paste0(figdir, "/matbiomass_sigma_rdevs.png"), dpi = 300, height = 6, width = 8, units = "in")

# II. Recruitment ----

# ASA-estimated age-3 abundance: 

# sigmaM appears to have a moderate affect on estimates but trends are the same.
# For larger sigmaMs: when survival is estimated to be low, recruitment is
# higher (and vice versa).
M + geom_line(aes(y = age3), size = 0.8) +
  labs(x = "\nYear", y = "ASA age-3 abundance (millions)\n", colour = "sigmaM")

ggsave(paste0(figdir, "/ASA_age3_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

# sigma_rdevs has no affect on estimates - very surprising to me! This is a
# good thing. Recruitment estimates are stable.
R + geom_line(aes(y = age3), size = 0.8) +
  labs(x = "\nYear", y = "ASA age-3 abundance (millions)\n", colour = "sigma_rdevs")

ggsave(paste0(figdir, "/ASA_age3_sigma_rdevs.png"), dpi = 300, height = 6, width = 8, units = "in")

# Ricker-estimated age-3 abundance:

# Very different trends in first part of ts, similar after 2000-ish. sigmaM
# affects both trend and magnitude of the estimates.
M + geom_line(aes(y = SR), size = 0.8) +
  labs(x = "\nYear", y = "Ricker age-3 abundance (millions)\n", colour = "sigmaM")

ggsave(paste0(figdir, "/Ricker_age3_sigmaM.png"), dpi = 300, height = 6, width = 8, units = "in")

# sigma_rdevs only mildly influential, and only at values of sigmaM=0.02 and
# 0.03, where higher sigma_rdevs result in lower estimated recruitment
  labs(x = "\nYear", y = "Ricker age-3 abundance (millions)\n", colour = "sigma_rdevs")

ggsave(paste0(figdir, "/Ricker_age3_sigma_rdevs.png"), dpi = 300, height = 6, width = 8, units = "in")
