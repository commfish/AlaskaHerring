
# Sensitivity analysis on sigmaM (variability block-specific natural mortality
# deviations) and conditioning on catch versus effort.
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov Last edited: 2019-07-09

# Description of code and variable definitions: Compare sensitivity results of
# sigmaM (the parameters controling variability of natural mortality deviations)
# with the 2018 best parameterization of LS and this model conditioned on
# effort.

# Pseudocode 
# 1. Set up directories (create a subdirectory of HER called sensitivity_sigmaM
# and copy the tpl to each run folder)
# 2. Select values of sigma_M over which to run sensitivity analysis
# 3. Use a nested for loop that fixes values of condition (on catch vs effort)
# over a range of sigmaM while keeping track of convergence and model results

# Set up/user inputs ----

YEAR <- 2019 # Forecast year

# IMPORTANT: User inputs model version name. This will create a subdirectory and
# allow you to run multiple model versions for comparison.
MODEL_VERSION <- "HER_bestLS_321"   # HER with the "best" LS parameterization

source(paste0(YEAR, "_forecast/r/create_ctl.r")) # function that writes new ctl file
source(paste0(YEAR, "_forecast/r/helper.r"))

options(scipen = 999) # remove scientific notation

root_dir <- getwd() # project root
tpl_dir <- file.path(root_dir, paste0(YEAR, "_forecast/admb/", MODEL_VERSION)) # location of tpl
run_dir <- file.path(root_dir, paste0(YEAR, "_forecast/results/sensitivity_sigmaM"))
dir.create(run_dir, showWarnings = FALSE)
files <- list.files(file.path(tpl_dir)) # searchable file list
tpl <- file.path(tpl_dir, files[which(grepl(".tpl", files))]) # find tpl
name <- strsplit(files[grepl("tpl",files)],".tpl")[[1]][1] # program name
dat <- file.path(tpl_dir, files[which(grepl(".dat", files))]) # find dat(two)

# Compile tpl if it hasn't been compiled yet
setwd(tpl_dir)
exe <- paste0(name,".exe")
if (file.exists(exe) == FALSE) {
  setup_admb() # check admb path
  compile_admb(name, verbose = TRUE) # compile
} 

std <- paste0(name,".std") # for checking model convergence in loop

# Run analysis ----

# Variability of recruitment devs defined by precision = precR = 1/(sigmaR^2) or
# sigmaR = sqrt(1/precR). Smaller precisions result in large sigmaR and vice
# versa. Results from 2018_forecast/text/Sensitivity on sigmaR, sigmaM, and
# sigma_rdevs.doc recommend sigmaR=0.5. Fix sigmaR at this value.
sigmaR <- 0.5
precR <- round(1/(sigmaR^2), 3)

# Condition - on catch or effort? catch = 0, effort = 1
condition_vec <- c(0, 1)
condition_name <- c("Catch", "Effort")

# Results from 2018_forecast/text/Sensitivity on sigmaR, sigmaM, and
# sigma_rdevs.doc recommend searching over 0.05-0.08. Slightly expand range
# here:
sigmaM_vec <- seq(0.05, 0.1, 0.01) 

# check_converge <- list()
check_converge <- matrix(nrow = length(condition_vec), ncol = length(sigmaM_vec))

for(i in 1:length(condition_vec)){
  for(j in 1:length(sigmaM_vec)){

    sensdir <- file.path(run_dir, paste0("condition_", condition_name[i], "_sigmaM_", sigmaM_vec[j]))
    dir.create(sensdir, showWarnings = FALSE)
    setwd(sensdir)
    
    # Create new sitka_2.ctl file
    create_ctl(precR = precR, sigmaM = sigmaM_vec[j], condition = condition_vec[i])

    # copy TPL to run directory
    file.copy(from = tpl, to = sensdir, overwrite = TRUE)	
    
    # copy exe
    file.copy(from = file.path(tpl_dir, exe), to = sensdir, overwrite = TRUE)	
    
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
    
    # Also check max gradient component
    mgc <- readMGC(name) 
    if (mgc < 1e-3) {
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
  mutate(condition = condition_name) %>% 
  gather("sigmaM", "converged?", -condition) %>% 
  mutate(fn = apply(expand.grid(paste0("condition_", condition_name), 
                                paste0("sigmaM_", sigmaM_vec)), 
                    1, paste, collapse = "_"))  -> check

# Save relevant files for later use 
setwd(run_dir)

if(exists("check")) {
  write_csv(check, "sensitivity_analysis_results.csv")
  } else {
    check <- read_csv("sensitivity_analysis_results.csv") }

# 1. Characterize models that did not converge: All models conditioned on effort
# converged. Only 1/3 of models conditioned on catch converged (sigmaM = 0.06 and 0.09)
check %>% arrange(condition)

# Read results for converged models
filter(check, `converged?` == 1) -> check

diags_ls <- list() # store diagnostic output
derived_ls <- list() # store derived time series of interest
natmat_ls <- list() # time-varying natural mortality & survival
matsel_ls <- list() # time-varying maturity and selectivity

# NOTE- Year for output will some times indicate terminal year in a time block
for(i in 1:length(check$fn)){
  
  condition <- check$condition[i]
  sigmaM <- check$sigmaM[i]
  fn <- check$fn[i]
  
  sensdir <- file.path(run_dir, fn)
  setwd(sensdir)
  
  D <- read_admb(name)
  
  # Diagnostics, key derived forecast quantities and parameter estimates
  diags_ls[[i]] <- data.frame(fn = fn,
                              condition = condition,
                              sigmaM = sigmaM,
                              nll = D[["fit"]]$nlogl,
                              maxgrad = D[["fit"]]$maxgrad,
                              nopar = D[["fit"]]$nopar,
                              ghl = D[["ghl"]], # Already in short tons
                              fore_sb = D[["fore_sb"]] / 0.90718, # convert to short tons
                              fore_matb = D[["fore_matb"]] / 0.90718, # convert to short tons
                              Mbar = exp(D[["theta"]][1]),
                              rinit = exp(D[["theta"]][2]),
                              rbar = exp(D[["theta"]][3]))
  
  # Natural mortality and survival
  df <- data.frame(D[["year"]], 
                   D[["Mij"]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  
  natmat_ls[[i]] <- df %>% 
    select(Year, M = `3`) %>% # nat mat not age-specific, so just use one age
    mutate(survival = exp(-M),
           condition = condition,
           sigmaM = sigmaM,
           fn = fn) %>% 
    filter(Year %in% c(1998, 2014, 2018)) %>% # Predefined time blocks *FLAG* this is hard coded
    mutate(`Time blocks` = c("1980-1998", "1999-2014", "2015-2018"),
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
                        labels = c("3", "4", "5", "6", "7", "8+"))) %>% 
    group_by(Age, param, proportion) %>% 
    mutate(min = min(Year),
           max = max(Year),
           `Time blocks` = paste0(min, "-", max),
           condition = condition,
           sigmaM = sigmaM,
           fn = fn) %>%
    ungroup() %>% 
    distinct(condition, sigmaM, fn, param, `Time blocks`, Year = max, Age, proportion)
  
  # Derived time series of biomass and ASA vs Ricker-predicted recruitment
  derived_ls[[i]] <- data.frame(Year = D[["year"]], 
                                spB = D[["sp_B"]] / 0.90718, # convert to short tons
                                matB = D[["mat_B"]] / 0.90718,
                                age3 = D[["Nij"]][1:nyr, 1] # ASA age-3 recruits
  ) %>% 
    left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                    D[["mod_nyr"]] + 1, 1),
                         # SR is the current year age-3
                         SR = D[["recruits"]],
                         # ssb is the spawning biomass in year-3 that produced
                         # the current years recruits
                         ssb = D[["spawners"]] / 0.90718), by = "Year") %>% 
    mutate(condition = condition,
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

# For showing conditioned of catch vs effort
palR <- colorRampPalette(c("#bfd3e6","#8c6bb1", "#4d004b"))
condition_cols <- palR(length(condition_vec))

# For showing sigmaM gradient
pal <- colorRampPalette(c("#edf8b1","#7fcdbb", "#2c7fb8"))
sigmaM_cols <- pal(length(sigmaM_vec))

# Diagnostics ----

setwd(run_dir)
figdir <- file.path(run_dir, "figures")
dir.create(figdir, showWarnings = FALSE)

# Maximum gradient component
diags %>% 
  ggplot(aes(x = sigmaM, y = condition, size = maxgrad)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Maximum\ngradient",
             breaks = seq(min(diags$maxgrad), max(diags$maxgrad), length.out = 3),
             labels = round(seq(min(diags$maxgrad), max(diags$maxgrad), length.out = 3), 5)) +
  geom_point(data = filter(diags, maxgrad == min(maxgrad)),
             col = "green", size = 10, fill = NA, shape = 6) +
  geom_point(data = filter(diags, maxgrad == max(maxgrad)),
             col = "red", size = 10, fill = NA, shape = 2) +
  labs(x = "\nsigmaM", y = "Condition\n")

ggsave(paste0(figdir, "/maxgrad.png"), dpi = 300, height = 4, width = 5, units = "in")

# Objective function - nothing about this is surprisingly or even that
# informative. If estimated, these parameters would shrink to 0 because less
# variability is going to result in a better fit. The fact that this doesn't
# result in a huge spread of NLL values is a good sign.
diags %>% 
  ggplot(aes(x = sigmaM, y = condition, size = nll)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(trans = "reverse",
             range = c(0, 4),
             name = "Negative\nlog\nlikelihood",
             breaks = seq(min(diags$nll), max(diags$nll), length.out = 3),
             labels = round(seq(min(diags$nll), max(diags$nll), length.out = 3), 0)) +
  labs(x = "\nsigmaM", y = "Condition\n")

ggsave(paste0(figdir, "/nll.png"), dpi = 300, height = 4, width = 5, units = "in")

# Mean natural mortality (Mbar) - Very little contrast
diags %>% 
  ggplot(aes(x = sigmaM, y = condition, size = Mbar)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Mean natural\nmortality",
             breaks = seq(min(diags$Mbar), max(diags$Mbar), length.out = 5),
             labels = round(seq(min(diags$Mbar), max(diags$Mbar), length.out = 5), 3)) +
  labs(x = "\nsigmaM", y = "Condition\n")

ggsave(paste0(figdir, "/Mbar.png"), dpi = 300, height = 4, width = 5, units = "in")


# Initial mean numbers-at-age (rinits) - There is very little contrast in rinits
diags %>% 
  ggplot(aes(x = sigmaM, y = condition, size = rinit)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Initial mean\nnumbers-at-age\n(millions)",
             breaks = seq(min(diags$rinit), max(diags$rinit), length.out = 5),
             labels = round(seq(min(diags$rinit), max(diags$rinit), length.out = 5), 0)) +
  labs(x = "\nsigmaM", y = "Condition\n")

ggsave(paste0(figdir, "/rinit.png"), dpi = 300, height = 4, width = 5, units = "in")

# Age-3 mean recruitment (rbar) - Exact same results and rationale as initial
# mean numbers-at-age
diags %>% 
  ggplot(aes(x = sigmaM, y = condition, size = rbar)) +
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4),
             name = "Age-3\nmean\nrecruitment\n(millions)",
             breaks = seq(min(diags$rbar), max(diags$rbar), length.out = 5),
             labels = round(seq(min(diags$rbar), max(diags$rbar), length.out = 5), 0)) +
  labs(x = "\nsigmaM", y = "Condition\n")

ggsave(paste0(figdir, "/rbar.png"), dpi = 300, height = 4, width = 5, units = "in")

# GHL ----

# Large sigmaM's result in drastically smaller GHLs overall. The minimum,
# maximum, and approximate median values are shown in red triangle, green
# upsidedown triangle, and orange square.

diags %>% 
  ggplot(aes(x = sigmaM, y = condition, size = ghl)) +
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
  labs(x = "\nsigmaM", y = "Condition\n")

ggsave(paste0(figdir, "/ghl.png"), dpi = 300, height = 4, width = 5, units = "in")

# Natural mortality ----

# increasing sigmaM has a relatively small effect on time-varying M over this
# range
ggplot(natmat, aes(x = `Time blocks`, y = survival, 
                   group = sigmaM, colour = sigmaM)) +
  geom_line(size = 1) + 
  geom_point(size = 1.3) + 
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~condition) +
  labs(colour = "sigmaM", y = "Survival")

ggsave(paste0(figdir, "/survival_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

# For the two sigmaM values that converged when model was conditioned on catch,
# the survival estimates were the same regardless of conidition.
check %>% filter(condition == "Catch") %>% pull(sigmaM) -> sigmaM_compare

ggplot(natmat %>% 
         filter(sigmaM %in% sigmaM_compare), aes(x = `Time blocks`, y = survival, 
                   group = condition, colour = condition,
                   linetype = condition)) +
  geom_line(size = 1) + 
  geom_point(size = 1.3) + 
  scale_colour_manual(values = condition_cols) +
  facet_wrap(~sigmaM) +
  labs(colour = "Condition", linetype = "Condition", y = "Survival")

ggsave(paste0(figdir, "/survival_condition.png"), dpi = 300, height = 4, width = 8, units = "in")

# Maturity and selectivity ----
matsel %>% 
  group_by(param) %>% 
  mutate(combos = paste0("sigmaM: ", sigmaM, " (", `Time blocks`, ")"),
         combos_cond = paste0("Conditioned on ", condition, " (", `Time blocks`, ")")) -> matsel

# Trends appear constant across condition Increasing sigmaM appears to result in
# the maturity curves shifting slightly to the left (maturing earlier)
matsel %>% filter(param == "Maturity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = sigmaM, 
                group = combos), size = 1) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~condition) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Maturity\n", colour = "sigmaM")

ggsave(paste0(figdir, "/maturity_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

# Condition does not affect maturity estimates
matsel %>% filter(param == "Maturity" & sigmaM %in% sigmaM_compare) %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = condition, 
                group = combos_cond), size = 1) +
  scale_colour_manual(values = condition_cols) +
  facet_wrap(~sigmaM) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Maturity\n", colour = "Condition")

ggsave(paste0(figdir, "/maturity_condition.png"), dpi = 300, height = 4, width = 8, units = "in")

# sigmaM and condition do not influence selectivity estimates
matsel %>% filter(param == "Selectivity") %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, colour = sigmaM, 
                group = combos), size = 1) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~condition) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Selecitivity\n", colour = "sigmaM")

ggsave(paste0(figdir, "/selectivity_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

matsel %>% filter(param == "Selectivity" & sigmaM %in% sigmaM_compare) %>% 
  ggplot(aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = condition, colour = condition, 
                group = combos_cond), size = 1) +
  scale_colour_manual(values = condition_cols) +
  facet_wrap(~sigmaM) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  labs(x = "\nAge", y = "Selecitivity\n", colour = "Condition", linetype = "Condition")

ggsave(paste0(figdir, "/selectivity_condition.png"), dpi = 300, height = 4, width = 8, units = "in")

# Derived time series ----

# Separate forecast quantities into separate df - cannot use the model-supplied
# fore_sb b/c the assumptions are not correct
diags %>% 
  select(fn, condition, sigmaM, fore_matb, ghl) %>% 
  mutate(Year = max(derived$Year) + 1,
         fore_sb = fore_matb - ghl) -> forec

# SSB facetted by sigmaM - shows variability in ssb is attributed to sigmaM not
# condition

# Facet by condition (colours are sigmaM)
M <- derived %>% 
  ggplot(aes(x = Year, colour = sigmaM, group = sigmaM)) +
  scale_colour_manual(values = sigmaM_cols) +
  facet_wrap(~condition) +
  scale_y_continuous(label = scales::comma)

# Facet by sigmaM (colours are condition)
R <- derived %>% filter(sigmaM %in% sigmaM_compare) %>% 
  ggplot(aes(x = Year, colour = condition, linetype = condition, group = condition)) +
  scale_colour_manual(values = condition_cols) +
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

# condition has no influence on estimates
R + geom_line(aes(y = spB), size = 0.8) +
  labs(x = "\nYear", y = "Spawning biomass (t)\n", colour = "Condition", linetype = "Condition") +
  geom_point(data = forec %>% filter(sigmaM %in% sigmaM_compare), aes(x = Year, y = fore_sb, colour = condition, group = condition),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/spbiomass_condition.png"), dpi = 300, height = 4, width = 8, units = "in")

# Mature biomass: same results as spawning biomass
M + geom_line(aes(y = matB), size = 0.8) +
  labs(x = "\nYear", y = "Mature biomass (t)\n", colour = "sigmaM") +
  geom_point(data = forec, aes(x = Year, y = fore_matb, colour = sigmaM, group = sigmaM),
             shape = "*", size = 4)

ggsave(paste0(figdir, "/matbiomass_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

R + geom_line(aes(y = matB), size = 0.8) +
  labs(x = "\nYear", y = "Mature biomass (t)\n", colour = "Condition", linetype = "Condition") +
  geom_point(data = forec %>% filter(sigmaM %in% sigmaM_compare), aes(x = Year, y = fore_matb, colour = condition, group = condition),
             shape = "*", size = 4) 

ggsave(paste0(figdir, "/matbiomass_condition.png"), dpi = 300, height = 4, width = 8, units = "in")

# Recruitment ----

# ASA-estimated age-3 abundance: 

# sigmaM appears to have a slight affect on estimates but trends are the same.
M + geom_line(aes(y = age3), size = 0.8) +
  labs(x = "\nYear", y = "ASA age-3 abundance (millions)\n", colour = "sigmaM")

ggsave(paste0(figdir, "/ASA_age3_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

# condition has no affect on estimates
R + geom_line(aes(y = age3), size = 0.8) +
  labs(x = "\nYear", y = "ASA age-3 abundance (millions)\n", colour = "Condition", linetype = "Condition")

ggsave(paste0(figdir, "/ASA_age3_condition.png"), dpi = 300, height = 4, width = 8, units = "in")

# Ricker-estimated age-3 abundance:

# Very different trends in first part of ts, similar after 2000-ish.
M + geom_line(aes(y = SR), size = 0.8) +
  labs(x = "\nYear", y = "Ricker age-3 abundance (millions)\n", colour = "sigmaM")

ggsave(paste0(figdir, "/Ricker_age3_sigmaM.png"), dpi = 300, height = 4, width = 8, units = "in")

# Condition has small effect, especially in first part of time series. Trends
# over time are the same for a given sigmaM
R + geom_line(aes(y = SR), size = 0.8) +
  labs(x = "\nYear", y = "Ricker age-3 abundance (millions)\n", colour = "Condition", linetype = "Condition")

ggsave(paste0(figdir, "/Ricker_age3_condition.png"), dpi = 300, height = 4, width = 8, units = "in")

