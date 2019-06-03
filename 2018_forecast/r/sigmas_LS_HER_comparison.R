# Sensitivity of sigmaR and sigmaM with comparison figures to LS
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-06-03

# Compare sensitivity results of sigmaR and sigmaM (the parameters controling
# the variability of recruitment and time-varying natural mortality) with the
# 2017 best parameterization of LS. For both sigmaR and sigmaM compare low,
# medium, and high values (9 trials). Only create comparison figures for model
# runs that converged.

# Set up ----

# Forecast year
YEAR <- 2018

# Create directory for HER/LS comparison figs, HER figs, and LS figs
root_dir <- getwd() # project root
tpl_dir <- file.path(root_dir, paste0(YEAR, "_forecast/admb/HER_bestLS")) # location of tpl
run_dir <- file.path(root_dir, paste0(YEAR, "_forecast/results/sigmas_LS_HER_comparison"))
dir.create(run_dir, showWarnings = FALSE)
files <- list.files(file.path(tpl_dir)) # searchable file list
tpl <- file.path(tpl_dir, files[which(grepl(".tpl", files))]) # find .tpl
ctl <- file.path(tpl_dir, files[which(grepl(".ctl", files))]) # find .ctl
name <- strsplit(files[grepl("tpl",files)],".tpl")[[1]][1] # program name
dat <- file.path(tpl_dir, files[which(grepl(".dat", files))]) # fine .dat files (two)
std <- paste0(name,".std") # used to check model convergence
exe <- paste0(name,".exe") # executible file so it doesn't need to be recompiled each time

source(paste0(YEAR, "_forecast/r/tools.r"))
source(paste0(YEAR, "_forecast/r/helper.r"))
source(paste0(YEAR, "_forecast/r/create_ctl.r")) # function that writes new ctl file

# LS results ----

# LS YEAR forecast results
LS_forec <- read_csv(paste0(YEAR, "_forecast/data/LS_", YEAR, "forec_results.csv"))
LS_byage <- read_csv(paste0(YEAR, "_forecast/data/LS_", YEAR, "forec_results_byage.csv"))
LS_byyear <- read_csv(paste0(YEAR, "_forecast/data/LS_", YEAR, "forec_results_byyear.csv"))

# Recruitment (age-3 abundance) and associated residuals in LS_byyear should be
# NA in the first three years, it is currently 0.
LS_byyear %>% 
  mutate(SR = ifelse(year %in% c(1980, 1981, 1982), NA, SR),
         res_SR = ifelse(year %in% c(1980, 1981, 1982), NA, res_SR)) -> LS_byyear

# LS results of year minus one, two, and three year estimates and forecast of
# pre-fishery mature biomass. "Tidied" version of the data in
# S:\Region1Shared-DCF\Research\Herring-Dive Fisheries\Herring\Year 2018
# Forecasts\Forecast Models\Sitka\best_model.57_321
# FINAL\ADMB\sitka_graphics.ctl
LS_yearminus <- read_csv(paste0(YEAR, "_forecast/data/yearminus_matbio.csv"))

# in SPAWN folder, numbers from excel spread sheet for spawn deposition.
LS_byyear$surv_est_spB <- c(35000,30000,29500,23500,38500,31000,25000,46000,58500,27000,23000,23500,
                            43351,37150,14941,34990,40827,28611,34942,44554,57988,58756,40366,55769,
                            69907,101305,66111,84501,247088,110946,126230,161904,62518,103267,48561, 
                            58183,77973,46919)

# survey estimated mature biomass (add catch to survey estimated spawning biomass)
LS_byyear %>% 
  mutate(surv_est_matbio = surv_est_spB + tcb) -> LS_byyear

# b/c the output wasn't easily accessible in the LS summary csv's, I rewrote the
# threshold and ghl calculation for visualization purposes. The harvest_thresh and thresh_denom are from
# the 2018 sitka forecast model best model folder in sitka.ctl 

harvest_thresh <- 25000 # Harvest Threshold (tons)
thresh_denom <- 20000 # Threshold calculation denominator (tons)

LS_byyear %>% 
  mutate(# Defined in regulation, if below there is no fishery
    Threshold = ifelse(year <= 1996, 7500,
                       ifelse (year <= 2009, 20000, 25000))) -> LS_byyear

# GHL calculation:
LS_forec %>% 
  summarize(tot_mat_B_tons = sum(for_mat_baa_tons)) %>% 
  mutate(year = max(LS_byyear$year) + 1,
         threshold = (2 + 8 * tot_mat_B_tons / thresh_denom) / 100,
         hr = ifelse(tot_mat_B_tons < harvest_thresh, 0, 
                     ifelse(threshold >= 0.2, 0.2, threshold)),
         GHL = ifelse(hr == 0.2, 0.2 * tot_mat_B_tons,
                      threshold * tot_mat_B_tons),
         # Defined in regulation
         Threshold = ifelse(year <= 1996, 7500,
                            ifelse (year <= 2009, 20000, 25000)),
         `Spawning biomass forecast` = tot_mat_B_tons - GHL) -> LS_ghl

# Run sensitivity analysis ----

# Compile tpl if it hasn't been compiled yet
setwd(tpl_dir)
exe <- paste0(name,".exe")
if (file.exists(exe) == FALSE) {
  setup_admb() # check admb path
  compile_admb(name, verbose = TRUE) # compile
} 

# Do sensistivity analysis on three levels of variability (low, med, and high)
# for coarse comparison. Parameterization of sigmaR as precision.
sigmaR_vec <- c(0.5, 1.0, 1.5)
precR_vec <- round(1/(sigmaR_vec^2), 3)
sigmaM_vec <- c(0.001, 0.075, 0.500)

check_converge <- matrix(nrow = length(precR_vec), ncol = length(sigmaM_vec))

for(i in 1:length(precR_vec)){
  for(j in 1:length(sigmaM_vec)){
    
    sensdir <- file.path(run_dir, paste0("precR_", precR_vec[i], "_sigmaM_", sigmaM_vec[j]))
    dir.create(sensdir, showWarnings = FALSE)
    setwd(sensdir)
    
    # Create new sitka_2.ctl file
    create_ctl(precR = precR_vec[i], sigmaM = sigmaM_vec[j])
    
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

# Check convergence ----

check <- as.data.frame(check_converge, row.names = FALSE)
names(check) <- sigmaM_vec
check %>% 
  mutate(precR = precR_vec) %>% 
  gather("sigmaM", "converged?", -precR) %>% 
  mutate(fn = paste0("precR_", precR, "_sigmaM_", sigmaM)) -> check

# Save relevant files for later use 
setwd(run_dir)

if(exists("check")) {
  write_csv(check, "converged.csv")
} else {
  check <- read_csv("converged.csv") }

filter(check, `converged?` == 0) %>% arrange(precR)

# Read results for converged models
filter(check, `converged?` == 1) -> check

for(i in 1:length(check$fn)){
  
  precR <- check$precR[i]
  sigmaM <- check$sigmaM[i]
  fn <- check$fn[i]
  
  sensdir <- file.path(run_dir, fn)
  setwd(sensdir)
  
  D <- read_admb(name)
  
  # Create fig directory
  figdir <- file.path(sensdir, "figures")
  dir.create(figdir, showWarnings = FALSE)
  
  # Spawning biomass (post-fishery) ----

  tot_yrs <- D[["dat_nyr"]] - D[["dat_syr"]] + 1
  
  df <- data.frame(Year = D[["year"]],
                   matB = D[["mat_B"]] / 0.90718,
                   spB = D[["sp_B"]] / 0.90718, # convert to short tons
                   catch = D[["data_catch"]][10:47, 2]#[nyr + 1, tot_yrs, 1), 2] # just the column of catch, already in short tons
  ) %>% 
    mutate(matB2 = spB + catch,
           Model = "HER") %>% 
    select(-catch) %>% 
    gather("Biomass", "tons", -c(Model, Year)) %>% 
    bind_rows(data.frame(Year = LS_byyear$year,
                         matB = LS_byyear$tot_mat_B_tons,
                         spB = LS_byyear$tot_sp_B_tons,
                         Model = "LS") %>% 
                gather("Biomass", "tons", -c(Model, Year)))
  
  srv_index <- data.frame(Year = LS_byyear$year,
                          srv_spB = LS_byyear$surv_est_spB,
                          Model = "Historical survey index")
  
  df %>% filter(Biomass %in% c("spB")) -> df
  
  axisx <- tickr(df, Year, 5)
  ggplot() +
    geom_line(data = df, aes(x = Year, y = tons, colour = Model, linetype = Model), size = 1) +
    geom_point(data = srv_index, aes(x = Year, y = srv_spB, shape = "Historical estimates from survey")) +
    scale_shape_manual(values = 1) +
    scale_colour_grey() +
    theme(legend.position = c(0.25, 0.7)) +
    scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "", y = "Spawning biomass (tons)\n", shape = "Data") -> spbiomass_plot
  
  ggsave(paste0(figdir, "/compare_spbiomass_plot.png"), plot = spbiomass_plot, dpi = 300, height = 4, width = 6, units = "in")
  
  # Recruitment ----
  
  # Comparison: 
  
  # Age-3 recruits compared recruitment estimated from Ricker
  df <- data.frame(Year = D[["years"]],
                   age3 = D[["Nij"]][, 1]) %>% 
    left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                    D[["mod_nyr"]] + 1, 1),
                         SR = D[["recruits"]],
                         ssb = D[["spawners"]],
                         resids = D[["resd_rec"]]), by = "Year") %>% 
    mutate(Model = "HER") %>% 
    bind_rows(
      LS_byyear %>% 
        select(Year = year, age3 = init_age_3, SR, resids = res_SR) %>% 
        mutate(Model = "LS"))
  
  axisr <- tickr(df, Year, 5)
  ggplot(df, aes(x = Year, colour = Model, shape = Model, linetype = Model)) +
    geom_line(aes(y = SR)) +
    geom_point(aes(y = age3), alpha = 0.9) +
    scale_colour_grey() +
    scale_x_continuous(breaks = axisr$breaks, labels = axisr$labels) +
    scale_y_continuous(labels = scales::comma) + 
    theme(legend.position = c(0.1, 0.85)) +
    labs(x = "\nYear", y = "Age-3 recruits (millions)\n") -> recruits
  
  ggsave(paste0(figdir, "/compare_recruit_plot.png"), plot = recruits, dpi = 300, height = 4, width = 6, units = "in")
  
  # Compare residuals
  ggplot(df, aes(x = Year, y = resids)) + 
    geom_hline(yintercept = 0, colour = "grey", size = 1) +
    geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
                 size = 0.2, colour = "grey") +
    geom_point() +
    facet_wrap(~ Model, ncol = 1) +
    labs(x = "\nYear", y = "Residuals\n") +
    scale_x_continuous(breaks = axisr$breaks, labels = axisr$labels) -> rec_resids
  
  ggsave(paste0(figdir, "/compare_recruitresids.png"), plot = rec_resids, dpi = 300, height = 5, width = 6, units = "in")
  
  # compare spawner recruit curves
  LS_byyear %>% 
    select(Year = year) %>% 
    filter(Year >= min(Year) + 3) %>% 
    bind_cols(LS_byyear %>% 
                # shift spawners up by 3 years 
                filter(year <= max(year) - 3) %>% 
                select(tot_sp_B_tons)) %>% 
    mutate(Model = "LS") %>% 
    right_join(df) %>% 
    mutate(ssb = ifelse(is.na(ssb), tot_sp_B_tons, ssb)) -> df
  
  ggplot(df, aes(x = ssb, colour = Model, shape = Model, linetype = Model)) +
    geom_line(aes(y = SR), size = 1) +
    geom_point(aes(y = age3)) +
    scale_colour_manual(values = c("black", "darkgrey")) + 
    theme(legend.position = c(0.1, 0.85)) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    geom_text_repel(aes(y = age3, label = Year), size = 2) +
    labs(x = "\nSpawning stock biomass (tons)", 
         y = "Age-3 recruits (millions)\n") -> sr_curve
  
  ggsave(paste0(figdir, "/compare_srcurves.png"), plot = sr_curve, dpi = 300, height = 6, width = 7, units = "in")
  
  # Egg deposition ----
  
  # Comparison:
  
  # From S. Dressel's bootstrap, have to add new values each year
  U <- c(1.18, 1.12, 1.10, 0.95, 1.23, 1.13, 0.93, 1.67, 2.51, 0.98, 0.80, 1.259686251, 1.851147636, 1.751126689, 0.560987576, 1.508697833, 1.633749193, 1.695401525, 1.255367509, 2.387620843, 2.795968238, 2.215761696, 1.462234716, 2.289501604, 2.650921062, 6.384923885, 2.279693146, 2.872760889, 29.05686308, 3.863145759, 4.816134565, 4.222205391, 1.634805164, 4.043867944, 1.288746439, 1.332721825, 2.445122264,1.69865191)
  L <- c(1.18, 1.12, 1.10, 0.95, 1.23, 1.13, 0.93, 1.67, 2.51, 0.98, 0.80, 0.985282715, 1.639191663, 1.382705136, 0.48874792, 1.398489625, 1.265029243, 1.286735024, 1.146877561, 1.827147032, 2.534746454, 1.882753246, 1.475607608, 1.863883108, 2.277982827, 3.540565615, 1.707320794, 2.568958439, 14.54490887, 3.237073047, 3.7630828, 3.942674884, 1.578639917, 2.996229014, 1.089460882, 1.154768444, 1.979792303, 1.357551379)
  
  df <- as.data.frame(D[["data_egg_dep"]][10:47, ])#[seq(tot_yrs - nyr + 1, tot_yrs, 1), ]) 
  colnames(df) <- c("Year", "obs", "log_se")
  df %>% filter(Year >= D[["mod_syr"]]) %>% 
    mutate(upper = obs + U,
           lower = obs - L) %>% 
    bind_cols(data.frame(HER = D[["pred_egg_dep"]],
                         LS = LS_byyear$tot_est_egg[1:38]))%>% 
    gather("Model", "trillions", -c(Year, obs, log_se, upper, lower)) -> df
  
  axisx <- tickr(df, Year, 5)
  ggplot(df, aes(x = Year)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), colour = "darkgrey", size = 0.001) + 
    geom_line(aes(y = trillions, colour = Model,
                  linetype = Model), size = 1) +  
    geom_point(aes(y = obs, shape = "Historical estimates from survey")) +
    scale_shape_manual(values = 1) +
    scale_color_grey() +
    theme(legend.position = c(0.25, 0.7)) +
    scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
    labs(x = "", y = "Eggs spawned (trillions)\n", shape = "Data") -> eggdep_plot
  
  ggsave(paste0(figdir, "/compare_eggdep_plot.png"), plot = eggdep_plot, dpi = 300, height = 4, width = 6, units = "in")
  
  # Compare resids
  data.frame(Year = D[["year"]],
             resids = D[["resd_egg_dep"]],
             Model = "HER") %>% 
    bind_rows(LS_byyear %>% 
                select(Year = year, resids = res_tot_egg) %>% 
                mutate(Model = "LS")) -> df
  
  axisx <- tickr(df, Year, 5)
  ggplot(df, aes(x = Year, y = resids)) + 
    geom_hline(yintercept = 0, colour = "grey", size = 1) +
    geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
                 size = 0.2, colour = "grey") +
    geom_point() +
    facet_wrap(~ Model, ncol = 1) +
    labs(x = "\nYear", y = "Residuals\n") +
    scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) -> egg_resids
  
  ggsave(paste0(figdir, "/compare_eggresids.png"), plot = egg_resids, dpi = 300, height = 5, width = 6, units = "in")
  
  # Comparison ----
  
  df <- data.frame(D[["year"]], 
                   D[["Mij"]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  write.table(df, sep="\t", row.names=FALSE)
  
  # Table summary
  df %>%
    select(Year, M = `3`) %>%
    mutate(HER = exp(-M)) %>% 
    select(-M) %>% 
    left_join(LS_byage %>% 
                select(Year, LS = survival), by = "Year") %>% 
    group_by(LS, HER) %>% 
    mutate(Blocks = paste0(min(Year), "-", max(Year))) %>% 
    distinct(Blocks, HER, LS) %>% 
    select(Blocks, HER, LS) -> survival_blks
  
  survival_blks[c(2:3)] <- lapply(survival_blks[,c(2:3)], prettyNum, digits = 4, big.mark=",")
  
  # Add in LS estimates
  df %>%
    select(Year, M = `3`) %>% 
    mutate(survival = exp(-M),
           Model = "HER") %>% 
    bind_rows(LS_byage %>% 
                select(Year, survival) %>% 
                mutate(Model = "LS")) -> df
  
  # Figure summary
  tickr(LS_byage, Year, 5) -> axis
  ggplot(df, aes(x = Year, y = survival,  colour = Model, linetype = Model)) +
    geom_vline(xintercept = c(1998.5, 2014.5), colour = "lightgrey", linetype = 3, alpha = 0.4) +
    geom_line(size = 1) +
    # geom_point() +
    lims(y = c(0, 1)) +
    annotation_custom(tableGrob(survival_blks, rows = NULL, 
                                theme = ttheme_minimal(base_size = 8, base_colour = "black", base_family = "Times",
                                                       parse = FALSE, #core = list(fg_params = list(hjust = 1)),
                                                       # colhead = list(fg_params = list(hjust = 1)), 
                                                       padding = unit(c(4, 4), "mm"))), 
                      xmin = 1980, xmax = 2015, ymin = 0.1, ymax = 0.35) +
    scale_colour_grey() +
    scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
    labs(x = "", y = "Survival\n") +
    theme(legend.position = c(0.1, 0.8)) -> survival_plot
  
  ggsave(paste0(figdir, "/compare_survival.png"), plot = survival_plot, dpi = 300, height = 4, width = 6, units = "in")
  
  # Maturity/Selectivity ----
  
  # Comparison:
  
  # First HER data
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
  
  df %>% 
    gather("Age", "proportion", -c(Year, param)) %>% 
    mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+")),
           # To make sure selectivity is differentiable, it was scaled to have a
           # mean of 1 across all ages. This was done in log space by substracting
           # the mean from the vector of age-specific selectivities. See Tech Doc
           # p 11. Here we normalize it from 0 to 1.
           proportion = ifelse(param == "Selectivity", (proportion - 0)/(max(proportion) - 0), 
                               proportion)) %>% 
    group_by(Age, param, proportion) %>% 
    mutate(min = min(Year),
           max = max(Year),
           `Time blocks` = paste0(min, "-", max)) %>% 
    group_by(param) %>% 
    mutate(Model = "HER",
           combos = paste0(Model, " ", `Time blocks`)) -> her_matsel
  
  # LS data 
  LS_byage %>% 
    select(Year, Age, Maturity = maturity, Selectivity = gear_select) %>% 
    mutate(Model = "LS") %>% 
    gather("param", "proportion", -c(Year, Age, Model)) %>% 
    group_by(Age, param, proportion) %>% 
    mutate(min = min(Year),
           max = max(Year),
           `Time blocks` = paste0(min, "-", max)) %>% 
    group_by(param) %>% 
    mutate(combos = paste0(Model, " ", `Time blocks`)) -> ls_matsel
  
  bind_rows(ls_matsel, her_matsel) -> matsel
  
  # In order to get separate legends for time block, this is split into separate
  # figures and then recombined.
  matsel %>% filter(param == "Maturity") -> par
  
  ggplot(par, aes(x = Age, y = proportion)) + 
    geom_line(aes(linetype = `Time blocks`, colour = Model, group = `combos`)) +
    geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
    expand_limits(y = 0) +
    scale_color_grey() +
    labs(x = "\nAge", y = "Proportion\n", linetype = "Time blocks") +
    ggtitle(paste0(par$param[1])) +
    theme(legend.position = c(0.7, 0.25),
          legend.spacing.y = unit(0, "cm"),
          plot.title = element_text(hjust = 0.5)) -> mat
  
  matsel %>% filter(param == "Selectivity") -> par
  
  ggplot(par, aes(x = Age, y = proportion)) + 
    geom_line(aes(linetype = `Time blocks`, colour = Model, group = `combos`)) +
    geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
    expand_limits(y = 0) +
    scale_color_grey() +
    labs(x = "\nAge", y = NULL, linetype = "Time blocks") +
    scale_y_continuous(breaks = seq(0, max(par$proportion), .25)) +
    ggtitle(paste0(par$param[1])) +
    theme(legend.position = c(0.7, 0.25),
          legend.spacing.y = unit(0, "cm"),
          plot.title = element_text(hjust = 0.5)) -> sel
  
  cowplot::plot_grid(mat, sel, align = "h", nrow = 1) -> matsel_plot
  
  ggsave(paste0(figdir, "/compare_mat_sel.png"), plot = matsel_plot, dpi = 300, height = 4, width = 6, units = "in")
  
  # Compare spawning comps ----
  
  df <- data.frame(spawnage_comp_obs = D[["data_sp_comp"]])
  colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
  df %>% 
    filter(Year >= 1980 ) %>%
    gather("Age", "proportion", -Year) %>% 
    mutate(Src = "Observed",
           Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+"))) -> obs
  
  df <- data.frame(Year = D[["year"]],
                   spawnage_comp_pred = D[["pred_sp_comp"]])
  colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
  df %>% 
    gather("Age", "proportion", -Year) %>% 
    mutate(Model = "HER",
           Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+"))) %>% 
    bind_rows(LS_byage %>% 
                select(Year, Age, proportion = spawnage_comp_est) %>% 
                mutate(Model = "LS") ) -> pred
  
  ggplot() +
    geom_bar(data = obs, aes(x = Age, y = proportion), 
             stat = "identity", colour = "grey", fill = "white", size = 0.4,
             width = 0.8, position = position_dodge(width = 0.5)) +
    geom_line(data = pred, aes(x = Age, y = proportion, colour = Model, 
                               linetype = Model, group = Model), size = 0.8) +
    facet_wrap(~ Year, dir = "v", ncol = 5) +
    # scale_colour_grey() +
    scale_colour_manual(values = c("black", "darkgrey")) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
    labs(x = '\nAge', y = 'Proportion-at-age\n') +
    theme(legend.position = "bottom") -> spcomp_barplot
  
  ggsave(paste0(figdir, "/compare_spcomp_barplot.png"), plot = spcomp_barplot, dpi = 300, height = 8, width = 6, units = "in")
  
  # Compare catch comps ----
  
  #commercial catch age comps, purse seine
  df <- data.frame(comage_comp_obs = D[["data_cm_comp"]])
  colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
  df %>% 
    filter(Year >= 1980 ) %>%
    gather("Age", "proportion", -Year) %>% 
    mutate(Src = "Observed",
           Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+"))) -> obs
  
  df <- data.frame(Year = D[["year"]],
                   comage_comp_pred = D[["pred_cm_comp"]])
  colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
  df %>% 
    gather("Age", "proportion", -Year) %>% 
    mutate(Model = "HER",
           Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                        labels = c("3", "4", "5", "6", "7", "8+"))) %>% 
    bind_rows(LS_byage %>% 
                select(Year, Age, proportion = catchage_comp_est) %>% 
                mutate(Model = "LS") ) -> pred
  
  ggplot() +
    geom_bar(data = obs, aes(x = Age, y = proportion), 
             stat = "identity", colour = "grey", fill = "white", size = 0.4,
             width = 0.8, position = position_dodge(width = 0.5)) +
    geom_line(data = pred, aes(x = Age, y = proportion, colour = Model, 
                               linetype = Model, group = Model), size = 0.8) +
    # scale_colour_grey() +
    scale_colour_manual(values = c("black", "darkgrey")) +
    facet_wrap(~ Year, dir = "v", ncol = 5) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
    labs(x = '\nAge', y = 'Proportion-at-age\n') +
    theme(legend.position = "bottom") -> catchcomp_barplot
  
  ggsave(paste0(figdir, "/compare_catchcomp_barplot.png"), plot = catchcomp_barplot, dpi = 300, height = 8, width = 6, units = "in")
  
}
