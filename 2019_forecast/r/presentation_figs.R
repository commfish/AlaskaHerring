# Run HER, including Bayesian analysis, and compare results to LS
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-07-15

# HER - ADFG's new herring model. Original ADMB code written by SJD Martell. Helper
# files and documentation contributed by M. Rudd.

# This code creates figures for the least-squares (LS) model and HER and a set
# of comparison figures. It runs the full Bayesian analysis.

# Set up and user inputs ----

set.seed(907) # random seed to make prediction interval calculations repeatable

library(ggdistribute) # plot posteriors
library(coda) # mcmc diagnostics
library(BayesianTools) # more mcmc diagnostics
# devtools::install_github('cttobin/ggthemr')
library(ggthemr)

# IMPORTANT: User inputs model version name. This will create a subdirectory and
# allow you to run multiple model versions for comparison.
# MODEL_VERSION <- "HER_bestLS_321" 
MODEL_VERSION <- "HER_best_condEffort.12_322"   # HER with best HER parameterization by AIC, conditioned on effort
# MODEL_VERSION <- "HER_best_condCatch.12_322"   # HER with best HER parameterization by AIC, conditioned on catch


# Forecast year
YEAR <- 2019

# User input: do you want to run the full MCMC?
run_mcmc <- FALSE

# Do you want to run LS figures? Option available so you don't end up with a
# bunch of repeat figures
run_LSfig <- FALSE

# Create directory for HER/LS comparison figs, HER figs, and LS figs
main_dir <- getwd()
# Directory for all comparison figures
fig_dir <- file.path(main_dir, paste0(YEAR, "_forecast/results/", MODEL_VERSION, "/presentation"))
dir.create(fig_dir, showWarnings = FALSE)

# source(paste0(YEAR, "_forecast/r/tools.r"))
source(paste0(YEAR, "_forecast/r/helper.r"))

my_cols <- c('#555555', '#008080', '#DAA520', '#330000') # first col is placeholder
ggthemr_reset() # remove previous effects
# Define colours for your figures with define_palette
tableau <- define_palette(
  swatch = my_cols, # colours for plotting points and bars
  gradient = c(lower = my_cols[1L], upper = my_cols[2L]), #upper and lower colours for continuous colours
  background =  "#f8f8f8" # "#F6F9FB"  #defining a grey-ish background 
)
# set the theme for your figures:
ggthemr(tableau, layout = "clean", type = "outer", text_size = 18)
# ggthemr(palette = "earth", spacing = 1, type = "outer", text_size = 20)

# LS results ----

# LS YEAR forecast results
LS_forec <- read_csv(paste0(YEAR, "_forecast/data/LS_", YEAR, "forec_results.csv"))
LS_byage <- read_csv(paste0(YEAR, "_forecast/data/LS_", YEAR, "forec_results_byage.csv"))
LS_byyear <- read_csv(paste0(YEAR, "_forecast/data/LS_", YEAR, "forec_results_byyear.csv"))

# Recruitment (age-3 abundance) and associated residuals in LS_byyear. Should be
# NA in the first three years, it is currently 0. 
LS_byyear %>% 
  mutate(SR = ifelse(SR == 0, NA, SR),
         res_SR = ifelse(res_SR == 0, NA, res_SR)) -> LS_byyear

LS_byyear %>% 
  select(-c(SR, res_SR)) %>% 
  full_join(LS_byyear %>% 
              select(year, SR, res_SR) %>% 
              mutate(year = year + 3) %>% 
              filter(year < max(year) - 2), by = "year") -> LS_byyear

# LS results of year minus one, two, and three year estimates and forecast of
# pre-fishery mature biomass. 
LS_yearminus <- read_csv(paste0(YEAR, "_forecast/data/yearminus_matbio.csv"))

# FLAG! The sitka and craig fig loop files have different definitions for 'num'
# (neither of which are accurate). numbers from excel spread sheet for spawn
# deposition in SPAWN folder. this 'spawn_dep' is NOT egg deposition, it's
# spawning biomass calculated in craig fig loop code as 1e12*tot_obs_egg/num,
# where num is not defined but is presumably spawners or weight or something
# like that.
LS_byyear$surv_est_spB <- LS_byyear$spawn_dep

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

# Run model in ADMB ----

# Can't figure out how to put file path directly from project root into the
# compile_admb() function
setwd(file.path(main_dir, paste0(YEAR, "_forecast/admb/", MODEL_VERSION)))

# Running model in R
setup_admb()
compile_admb("her", verbose = TRUE)
run_admb("her", verbose = TRUE) 

# Running her model in command line
# 1) Open command prompt ("C:/ADMB/ADMB Command Prompt (MinGW 64Bit)")
# 2) Navigate to AlaskaHerring deep inside S drive using 
# > cd ../..
# > S:
# > cd "S:\Region1Shared-DCF\Research\Herring-Dive Fisheries\Herring\ADMB Rudd Workshop 2018\AlaskaHerring\2018_forecast\admb\HER_bestLS"
# Compile
# > admb her 
# Run
# > her

# MCMC ----

# Number of interations (MCMC samples) and thinning rate (if 10, saves every
# 10th sample)
niter <- 11110000 # this value give you nout=10000
thin <- 1000
burn_in <- 0.1 # 10% (this is the first 10% of the thinned/saved iterations)
nout <- round((niter / thin + 1) - burn_in * (niter / thin + 1), 0)

if (run_mcmc == TRUE){
  run_admb("her", extra.args = paste0("-mcmc ", niter, " -mcsave ", thin))
  run_admb("her", extra.args="-mceval", verbose = TRUE)
}

D <- read_admb("her")

# NOTE: See helper.R for source code for all user-defined functions. All
# functions save a csv of the posterior/credible or posterior predictive
# interval in the HERfig_dir

# Posterior intervals (aka credible intervals, show variability around
# the expected value)
spb_sum <- ps_byyear(fn = "sp_B", unit_conv = 0.90718, save = FALSE)
matb_sum <- ps_byyear(fn = "mat_B", unit_conv = 0.90718, save = FALSE)
catch_sum <- ps_byyear(fn = "pred_catch", unit_conv = 0.90718, save = FALSE) # warning msg ok!
ricker_sum <- ps_byyear(fn = "ricker_rec", syr = D[["mod_syr"]] + D[["sage"]],
                         lyr = D[["mod_nyr"]] + 1, save = FALSE)
age3_sum <- ps_byyear(fn = "age3_rec", syr = D[["mod_syr"]],
                         lyr = D[["mod_nyr"]] + 1, save = FALSE)
fmort_sum <- ps_byyear(fn = "fmort", syr = D[["mod_syr"]], # won't be anything if conditioned on catch
                         lyr = D[["mod_nyr"]], save = FALSE)
egg_sum <- ps_byyear(fn = "pred_egg_dep", save = FALSE)

# Posterior predictive intervals for egg deposition
df <- as.data.frame(D[["data_egg_dep"]])
colnames(df) <- c("Year", "obs_egg", "log_se")
df %>% 
  filter(Year >= D[["mod_syr"]]) %>% 
  mutate(Year = factor(Year)) %>% 
  select(Year, log_se) %>% 
  left_join(egg_sum) %>%  
  mutate(log_value = log(value),
         # Posterior predictive value (draw from random lognormal distribution
         # given mean and sd of data)
         pp_value = exp(rnorm(n(), mean = log_value, sd = log_se))) %>% 
  group_by(Year) %>% 
  mutate(p025 = quantile(pp_value, 0.025),
         p975 = quantile(pp_value, 0.975),
         p250 = quantile(pp_value, 0.250),
         p750 = quantile(pp_value, 0.750)) %>% 
  distinct(Year, mean, median, p025, p975, p250, p750) %>% 
  ungroup() %>% 
  mutate(Year = as.numeric(as.character(Year))) -> egg_pp

# Survival posterior distribution
surv_sum <- ps_byage(fn = "survival", syr = D[["mod_syr"]],
                    lyr = D[["mod_nyr"]], n = niter / thin + 1, save = FALSE)

# Maturity posterior distribution
mat_sum <- ps_byage(fn = "maturity", syr = D[["mod_syr"]], 
                    lyr = D[["mod_nyr"]], n = niter / thin + 1, save = FALSE)

# Selectivity posterior distribution
sel_sum <- ps_byage(fn = "selectivity", syr = D[["mod_syr"]], 
                    lyr = D[["mod_nyr"]] + 1, n = niter / thin + 1, save = FALSE)  
sel_sum <- sel_sum[min != YEAR]

# Posterior intervals (credible intervals) for age compositions
sp_comp_sum <- ps_comps(fn = "pred_sp_comp", syr = D[["mod_syr"]], 
                        lyr = D[["mod_nyr"]], n = niter / thin + 1, save = FALSE) 

cm_comp_sum <- ps_comps(fn = "pred_cm_comp", syr = D[["mod_syr"]], 
                    lyr = D[["mod_nyr"]], n = niter / thin + 1, save = FALSE)  

# Posterior predictive intervals for age comp data

# Variance estimates for each posterior sample from the multivariate logistic
# distribution used for age comps
sp_tau2 <- ps_tau(fn = "pp_sp_tau2", save = FALSE)
cm_tau2 <- ps_tau(fn = "pp_cm_tau2", save = FALSE)

# Get posterior predictive intervals for multivariate logistic distribution (age
# compositions). Uses the posterior sample estimates and variance.
cm_comp_ppi <- ppi_rmvlogistic(fn = "cm_comp",
                               ps_sum = cm_comp_sum,
                               tau2 = cm_tau2, save = FALSE)

sp_comp_ppi <- ppi_rmvlogistic(fn = "sp_comp",
                               ps_sum = sp_comp_sum,
                               tau2 = sp_tau2, save = FALSE)

# Posterior distribution of forecasted mature biomass 
df <- data.frame(read.table("fore_matb.ps"))
colnames(df) <- "quantity"
df %>% 
  mutate(Year = D[["mod_nyr"]] + 1,
         quantity = quantity / 0.90718,
         iter = row_number(),
         mean = mean(quantity),
         median = median(quantity),
         # 95% cred
         q025 = quantile(quantity, 0.025),
         q975 = quantile(quantity, 0.975),
         # 50% cred
         q250 = quantile(quantity, 0.250),
         q750 = quantile(quantity, 0.750),
         # 5% cred
         q475 = quantile(quantity, 0.475),
         q525 = quantile(quantity, 0.525)) -> fore_matb_sum

# Mature biomass ----

tot_yrs <- D[["dat_nyr"]] - D[["dat_syr"]] + 1

# For HER with MCMC variance estimation:

# Combine estimates and forecast
matb_sum %>% 
  mutate(Year = as.numeric(as.character(Year)),
         type = "Estimate") %>% 
  bind_rows(fore_matb_sum %>% 
              mutate(type = "Forecast")) %>% 
  distinct(Year, mean, q025, q975, type) %>% 
  mutate(`Model estimates and forecast` = "HER with 95% credible interval") %>% 
  bind_rows(
    LS_byyear %>% 
      select(Year = year, mean = tot_mat_B_tons) %>% 
      mutate(type = "Estimate") %>% 
      bind_rows(
        LS_forec %>% 
          summarize(mean = sum(for_mat_baa_tons)) %>% 
          mutate(Year = max(LS_byyear$year) + 1,
                 type = "Forecast")) %>% 
      mutate(`Model estimates and forecast` = "LS")) %>% 
  mutate(label = ifelse(type == "Forecast", 
                        prettyNum(mean, big.mark = ",", digits = 1), NA)) -> df

tickryr <- data.frame(Year = 1980:2023)
axisf <- tickr(tickryr, Year, 5)
ggplot(df, aes(x = Year)) +
  geom_ribbon(data = filter(df, !is.na(q025)), aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.2, fill = '#008080') +
  geom_line(data = df %>% filter(type == "Estimate"), aes(y = mean, colour = `Model estimates and forecast`, linetype = `Model estimates and forecast`, group = `Model estimates and forecast`), size = 1) +
  geom_point(data = LS_byyear, aes(x = year, y = surv_est_matbio, shape = "Historical estimates from survey"), 
             colour = "black", size = 2) +
  scale_shape_manual(values = 1) +
  geom_point(data = df %>% filter(type == "Forecast"), aes(x = Year, y = mean, colour = `Model estimates and forecast`), shape = 8, size = 2) +
  geom_text_repel(data = df, aes(x = Year, y = mean, label = label, colour = `Model estimates and forecast`),
                  size = 5, nudge_x = 1.5, show.legend = FALSE) +
  scale_x_continuous(limits = c(1980, 2022), labels = seq(1980, 2020, 5), breaks = seq(1980, 2020, 5)) +
  scale_y_continuous(limits = c(0, max(LS_byyear$surv_est_matbio)), 
                     breaks = seq(0, max(LS_byyear$surv_est_matbio), 50000),
                     labels = scales::comma) +
  labs(x = NULL, y = NULL, #y = "Mature\nbiomass\n(tons)", 
       shape = NULL, text = NULL) +
  theme(legend.position = c(0.22, 0.7),
        axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  ggtitle("Mature biomass (tons)") -> post_matb

ggsave(paste0(fig_dir, "/mature_biomass.png"), plot = post_matb, dpi = 300, height = 5, width = 10, units = "in")

# Spawning biomass ----

# post-fishery
spb_sum %>% 
  mutate(Year = as.numeric(as.character(Year)),
         type = "Estimate") %>% 
  distinct(Year, mean, q025, q975, type) %>% 
  mutate(`Model estimates` = "HER with 95% credible interval") %>% 
  bind_rows(
    LS_byyear %>% 
      select(Year = year, mean = tot_sp_B_tons) %>% 
      mutate(type = "Estimate",
             `Model estimates` = "LS")) -> df

srv_index <- data.frame(Year = LS_byyear$year,
                        srv_spB = LS_byyear$surv_est_spB,
                        Model = "Historical estimates from survey")

ggplot() +
  geom_ribbon(data = filter(df, !is.na(q025)), aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.2, fill = '#008080') +
  geom_line(data = df, aes(x = Year, y = mean, colour = `Model estimates`, linetype = `Model estimates`), size = 1) +
  geom_point(data = srv_index, aes(x = Year, y = srv_spB, shape = "Historical estimates from survey"),
             colour = "black", size = 2) +
  scale_shape_manual(values = 1) +
  theme(legend.position = c(0.22, 0.7),
        axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  scale_y_continuous(limits = c(0, max(LS_byyear$surv_est_spB)), breaks = seq(0, max(LS_byyear$surv_est_spB), 50000), labels = scales::comma) +
  labs(x = "", y = NULL, #y = "Spawning\nbiomass\n(tons)", 
       shape = NULL) +
  ggtitle("Spawning biomass (tons)") -> spbiomass_plot

ggsave(paste0(fig_dir, "/spawning_biomass.png"), plot = spbiomass_plot, dpi = 300, height = 5, width = 10, units = "in")

# Abundance ----

tot_yrs <- D[["dat_nyr"]] - D[["dat_syr"]] + 1 # includes data
mod_yrs <- D[["mod_nyr"]] - D[["mod_syr"]] + 1 # just model years

# Spawning abundance: Convert sp_Nij to sp_Ni in millions
df <- D[["mat_Nij"]][1:mod_yrs,] # exclude the forecast year
df <- as.data.frame(df) 
colnames(df) <- c(paste(D[['sage']]:D[['nage']]))

mat_N <- df %>% mutate(year = (D[["mod_syr"]]):(D[["mod_nyr"]])) %>% 
  gather("Age", "mat_N", -year) %>% 
  group_by(year) %>% 
  summarize(mat_N = sum(mat_N))

# Total abundance: Convert Nij to Ni in millions
df <- D[["Nij"]][1:mod_yrs,] # exclude the forecast year
df <- as.data.frame(df) 
colnames(df) <- c(paste(D[['sage']]:D[['nage']]))

N <- df %>% mutate(year = (D[["mod_syr"]]):(D[["mod_nyr"]])) %>% 
  gather("Age", "N", -year) %>% 
  group_by(year) %>% 
  summarize(N = sum(N))

# Join and gather indices
mat_N %>% left_join(N) %>% 
  gather("var", "value", -c(year)) %>% 
  mutate(Model = "HER") %>% 
  bind_rows(LS_byyear %>%
              select(year, N, mat_N = tot_mat_N) %>% 
              gather("var", "value", -c(year)) %>% 
              mutate(Model = "LS")) -> df

# Mature and total abundance
ggplot(df %>% 
         filter(var %in% c("mat_N", "N")) %>% 
         mutate(var = factor(var, ordered = TRUE, 
                             levels = c("N", "mat_N"),
                             labels = c("Total", "Mature"))),
       aes(x = year)) +
  geom_line(aes(y = value, col = Model, linetype = var), size = 1) +
  scale_linetype_manual(values = c(1, 4)) +
  scale_colour_manual(values = c('#008080', '#DAA520')) +
  theme(legend.position = c(0.1, 0.75),
        axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  scale_y_continuous(labels = scales::comma, 
                     limits = c(0, max(N$N) * 1.4)) +
  labs(x = NULL, y = NULL,#y = "Age-3+\nabundance\n(millions)", 
       linetype = NULL, fill = NULL) +
  guides(linetype  = guide_legend(order = 2),
         color = guide_legend(order = 1)) +
  ggtitle("Age-3+ abundance (millions)") ->  n

ggsave(paste0(fig_dir, "/abundance.png"), plot = n, dpi = 300, height = 5, width = 10, units = "in")

# Fishing mortality ----

df <- data.frame(Year = D[["mod_syr"]]:D[["mod_nyr"]],
                 ft = D[["ft"]])
df %>% 
  mutate(Model = "HER") %>% 
  bind_cols(fmort_sum %>% distinct(Year, q025, q975)) -> df

ggplot(df, aes(x = Year)) +
  geom_ribbon(aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.2, fill = '#008080') +
  geom_line(aes(y = ft, colour = Model,
                linetype = Model), size = 1) +  
  geom_point(aes(y = ft, colour = Model), size = 2) +
  theme(legend.position = c(0.1, 0.85)) +
  guides(shape  = guide_legend(order = 2),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 1)) +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(x = "", y = NULL, shape = NULL) +
  ggtitle("Fishing mortality") -> fmort_plot
fmort_plot
ggsave(paste0(fig_dir, "/fmort.png"), plot = fmort_plot, dpi = 300, height = 5, width = 10, units = "in")

# Catch  ----

df <- as.data.frame(D[["data_catch"]])
colnames(df) <- c("Year", "obs", "log_se")
df %>% filter(Year >= D[["mod_syr"]]) %>% 
  mutate(obs = obs / 0.90718,
         ln_catch = log(obs),
         std = 1.96 * sqrt(log(log_se + 1)),
         upper = exp(ln_catch + std),
         lower = exp(ln_catch - std),
         pred = D[["pred_catch"]] / 0.90718,
         Model = "HER") %>% 
  bind_cols(catch_sum %>% distinct(Year, q025, q975)) -> df

ggplot(df, aes(x = Year)) +
  # geom_errorbar(aes(ymin = lower, ymax = upper),
  #               colour = "grey15", size = 0.1, width = 0.3) +
  geom_ribbon(aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.2, fill = '#008080') +
  geom_line(aes(y = pred, colour = Model,
                linetype = Model), size = 1) +  
  geom_point(aes(y = obs, shape = "Data"),
             colour = "grey15", size = 2) +
  scale_shape_manual(values = 1) +
  theme(legend.position = c(0.1, 0.8)) +
  guides(shape  = guide_legend(order = 2),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 1)) +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  scale_y_continuous(limits = c(0, max(df$q975)), labels = scales::comma) +
  labs(x = "", y = NULL, shape = NULL) +
  ggtitle("Catch (tons)") -> catch_plot
catch_plot
ggsave(paste0(fig_dir, "/catch.png"), plot = catch_plot, dpi = 300, height = 5, width = 10, units = "in")

# Compare resids
data.frame(Year = D[["year"]],
           resids = D[["resd_catch"]],
           Model = "HER") -> df

ggplot(df, aes(x = Year, y = resids, colour = Model, shape = Model)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point(size = 2) +
  # facet_wrap(~ Model, ncol = 1) +
  labs(x = NULL, y = NULL) +
  ggtitle("Catch residuals") +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), 
                     breaks = seq(1980, YEAR-1, 5)) +
  theme(legend.position = 'none') -> catch_resids
catch_resids
ggsave(paste0(fig_dir, "/catch_resids.png"), plot = catch_resids, dpi = 300, height = 5, width = 10, units = "in")

# Catch / spawning biomass barplot ----
LS_byyear %>%
  select(year, Catch = tcb, `Spawning biomass` = tot_sp_B_tons, Threshold) %>% 
  full_join(LS_ghl %>% select(year, `Spawning biomass forecast`, GHL, Threshold)) %>% 
  gather("var", "value", -c(year, Threshold)) %>% 
  mutate(var = factor(var, ordered = TRUE, 
                      levels = c("Spawning biomass", "Spawning biomass forecast", "Catch", "GHL"))) -> df
axisb <- tickr(df, year, 5)
ggplot(df, aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  geom_line(aes(y = Threshold, linetype = "Threshold"), size = 1, colour = "black") +
  scale_fill_manual(values = c("white", "grey85", "grey50", "grey20")) +
  theme(legend.position = c(0.2, 0.8),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 125000, 10000)) +
  labs(x = "", y = "Biomass (tons)\n", linetype = NULL, fill = NULL) -> catch_plot

df <- data.frame(year = D[["mod_syr"]]:D[["mod_nyr"]],
                 catch = D[["data_catch"]][10:48,2] / 0.90718,
                mat_B = D[["mat_B"]] / 0.90718,
                Threshold = LS_byyear$Threshold) %>% 
  mutate(sp_B = mat_B - catch) %>% 
  select(-mat_B)

fore_matb_sum %>% distinct(Year, mean, median, q025, q975, q250, q750) -> f_mat

forec <- data.frame(year = YEAR,
                    ghl = D[["ghl"]],
                    Threshold = 25000) %>% 
  mutate(sp_B_forec = f_mat$mean - ghl)

df <- full_join(df, forec) %>% 
  gather("var", "value", -c(year, Threshold)) %>% 
  mutate(var = factor(var, ordered = TRUE, 
                      levels = c("sp_B", "sp_B_forec", "catch", "ghl"),
                      labels = c("Mature biomass", "Mature biomass forecast", "Catch", "GHL")))

axisb <- tickr(df, year, 5)
ggplot(df, aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  geom_line(aes(y = Threshold, linetype = "Threshold"), size = 1, colour = "black") +
  scale_fill_manual(values = c("white", "grey85", "grey50", "grey20")) +
  geom_errorbar(data = matb_sum %>% distinct(Year, q025, q975) %>% 
                  mutate(Year = as.numeric(as.character(Year))),
                aes(x = Year, ymin = q025, ymax = q975),
                colour = "black", size = 0.001, width = 0) +
  geom_errorbar(data = f_mat, aes(x = Year, ymin = q025, ymax = q975), colour = "black", size = 0.001, width = 0) + 
  theme(legend.position = c(0.2, 0.8),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "", y = "Biomass (tons)\n", linetype = NULL, fill = NULL) -> catch_plot

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
ggplot(df, aes(x = Year)) +
  geom_point(aes(y = age3, colour = Model, shape = Model, linetype = Model), size = 2) +
  geom_line(aes(y = age3, colour = Model, shape = Model, 
                linetype = Model), size = 1) +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  scale_y_continuous(labels = scales::comma) + 
  theme(legend.position = c(0.1, 0.85)) +
  labs(x = NULL, y = NULL) + 
  ggtitle("Age-3 recruits (millions)") -> recruits
recruits
ggsave(paste0(fig_dir, "/compare_age3_recruitment.png"), plot = recruits, dpi = 300, height = 5, width = 10, units = "in")

df %>% 
  mutate(Model = factor(Model,
                        levels = c("HER", "LS"),
                        labels = c("HER with 95% credible interval", "LS"))) %>% 
         ggplot(aes(x = Year)) +
  geom_point(aes(y = age3, colour = Model, shape = Model, linetype = Model), size = 2) +
  geom_ribbon(data = age3_sum %>% distinct(Year, q025, q975) %>%
                mutate(Year = as.numeric(as.character(Year))),
              aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.2, fill = '#008080') +
  geom_line(aes(y = age3, colour = Model, shape = Model, 
                linetype = Model), size = 1) +
  # scale_colour_grey() +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  scale_y_continuous(labels = scales::comma) + 
  theme(legend.position = c(0.8, 0.85)) +
  labs(x = NULL, y = NULL) + 
  ggtitle("Age-3 recruits (millions)") -> recruits
recruits
ggsave(paste0(fig_dir, "/compare_age3_recruitment_CI.png"), plot = recruits, dpi = 300, height = 5, width = 10, units = "in")
       
# Compare residuals
ggplot(df, aes(x = Year, y = resids, colour = Model, shape = Model)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point(size = 2) +
  facet_wrap(~ Model, ncol = 1) +
  labs(x = NULL, y = NULL) +
  ggtitle("Recruitment residuals") +
  scale_x_continuous(limits = c(1980, YEAR-1), 
                     labels = seq(1980, YEAR-1, 5), 
                     breaks = seq(1980, YEAR-1, 5)) +
  theme(legend.position = 'none') -> rec_resids
rec_resids
ggsave(paste0(fig_dir, "/rec_resids.png"), plot = rec_resids, dpi = 300, height = 5, width = 10, units = "in")

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
  # scale_colour_manual(values = c("black", "darkgrey")) + 
  theme(legend.position = c(0.1, 0.85)) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  geom_text_repel(aes(y = age3, label = Year), 
                  size = 2, show.legend = FALSE) +
  labs(x = "\nSpawning stock biomass (tons)", 
       y = "Age-3 recruits (millions)\n") +
  ggtitle("Ricker stock-recruitment") -> sr_curve
sr_curve
ggsave(paste0(fig_dir, "/compare_srcurves.png"), plot = sr_curve, dpi = 300, height = 7, width = 10, units = "in")

ricker_sum %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  distinct(Year, q025, q975) %>% 
  left_join(df) -> df

# Add uncertainty
df %>% 
  mutate(Model = factor(Model,
                        levels = c("HER", "LS"),
                        labels = c("HER with 95% credible interval", "LS"))) %>% 
  ggplot(aes(x = ssb)) +
  geom_line(aes(y = SR, colour = Model, shape = Model, linetype = Model), size = 1) +
  geom_point(aes(y = age3, colour = Model, shape = Model, linetype = Model)) +
  geom_ribbon(data = filter(df, Model == "HER"), aes(x = ssb, ymin = q025, ymax = q975),
              alpha = 0.3, fill = '#008080') +
  theme(legend.position = c(0.25, 0.85)) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  geom_text_repel(aes(y = age3, label = Year), 
                  size = 2, show.legend = FALSE) +
  labs(x = "\nSpawning stock biomass (tons)", 
       y = "Age-3 recruits (millions)\n") +
  ggtitle("Ricker stock-recruitment") -> sr_curve
sr_curve
ggsave(paste0(fig_dir, "/compare_srcurves_CI.png"), plot = sr_curve, dpi = 300, height = 7, width = 10, units = "in")


# Age-3 recruits compared recruitment estimated from Ricker
df <- data.frame(Year = D[["years"]],
                 age3 = D[["Nij"]][, 1]) %>% 
  left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                  D[["mod_nyr"]] + 1, 1),
                       SR = D[["recruits"]],
                       resids = D[["resd_rec"]]), by = "Year")

# Add in posterior intervals
age3_sum %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  distinct(Year, q025, q975) %>% 
  left_join(df) -> df2

ggplot(df2, aes(x = Year)) +
  geom_bar(aes(y = age3),
           stat = "identity", 
           alpha = 0.5,
           width = 0.8, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q025, ymax = q975), colour = "black", size = 0.001, width = 0) + 
  # If you don't want to show Ricker line, comment out:
  # geom_line(aes(y = SR, linetype = "Ricker model"), size = 1) +
  theme(legend.position = c(0.1, 0.75),
        legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.8, 'lines'),
        legend.spacing.y = unit(0, "cm")) +
  # scale_colour_manual(values = '#DAA520') +
  scale_x_continuous(limits = c(1979, YEAR+1), labels = seq(1980, YEAR+1, 5), breaks = seq(1980, YEAR+1, 5)) +
  scale_y_continuous(labels = scales::comma) + 
  theme(legend.position = c(0.15, 0.75)) +
  labs(x = NULL, y = NULL, colour = NULL,
       linetype = NULL) +
  ggtitle("Age-3 recruits (millions)") -> recruits

ggsave(paste0(fig_dir, "/age3_recruitment_barsCI.png"), plot = recruits, dpi = 300, height = 5, width = 10, units = "in")

# Compare both models
df2 %>% select(Year, age3) %>% 
  mutate(Model = "HER") %>% 
  bind_rows(data.frame(Year = D[["mod_syr"]]:(YEAR-1),
            age3 = LS_byyear$init_age_3,
            Model = "LS")) -> df2
df2 %>% 
  filter(Year < YEAR) %>% 
  ggplot(aes(x = Year)) +
  geom_bar(aes(y = age3, colour = Model, fill = Model),
           stat = "identity",
           alpha = 0.5,
           width = 0.7,
           position = position_dodge(width = 0.8)) +
  theme(legend.position = c(0.1, 0.8)) +
  scale_x_continuous(limits = c(1979, YEAR), labels = seq(1980, YEAR, 5), breaks = seq(1980, YEAR, 5)) +
  scale_y_continuous(labels = scales::comma) + 
  labs(x = NULL, y = NULL, #colour = NULL, fill = NULL,
       linetype = NULL) +
  ggtitle("Age-3 recruits (millions)") -> recruits
recruits
ggsave(paste0(fig_dir, "/compare_age3_recruitment_bars.png"), plot = recruits, dpi = 300, height = 5, width = 10, units = "in")

# Egg deposition ----

# Comparison:
df <- as.data.frame(D[["data_egg_dep"]][10:tot_yrs, ])#[seq(tot_yrs - nyr + 1, tot_yrs, 1), ]) 
colnames(df) <- c("Year", "obs", "log_se")
df %>% filter(Year >= D[["mod_syr"]]) %>% 
  # From S. Dressel's bootstrap, have to add new values each year
  # 2008 survey egg estimates were extreme high and variable, pmin() helps
  # reduce the scale of the upper confidence interval for better visualization:
  mutate(upper = pmin(obs + LS_byyear$egg_upper, max(obs) * 1.2),
         lower = obs - LS_byyear$egg_lower) %>% 
  bind_cols(data.frame(HER = D[["pred_egg_dep"]],
                       LS = LS_byyear$tot_est_egg))%>% 
  gather("Model", "trillions", -c(Year, obs, log_se, upper, lower)) -> df

ggplot(df, aes(x = Year)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey15", size = 0.1, width = 0.3) + 
  geom_line(aes(y = trillions, colour = Model,
                linetype = Model), size = 1) +  
  geom_point(aes(y = obs, shape = "Historical estimates from survey with 95% confidence intervals"),
             colour = "grey15", size = 2) +
  scale_shape_manual(values = 1) +
  theme(legend.position = c(0.35, 0.75)) +
  guides(shape  = guide_legend(order = 2),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 1)) +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  labs(x = "", y = NULL, shape = NULL) +
  ggtitle("Eggs spawned (trillions)") -> eggdep_plot
eggdep_plot
ggsave(paste0(fig_dir, "/compare_eggdep.png"), plot = eggdep_plot, dpi = 300, height = 5, width = 10, units = "in")

# With uncertainty

# add credible and posterior predictive intervals
egg_sum %>% 
  distinct(Year, q025, q975) %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  left_join(df) %>% 
  left_join(egg_pp %>% 
  distinct(Year, p025, p975) %>% 
  mutate(Year = as.numeric(as.character(Year)))) -> df

df %>% 
  mutate(Model = factor(Model,
                        levels = c("HER", "LS"),
                        labels = c("HER with 95% credible interval", "LS"))) %>% 
  ggplot(aes(x = Year)) +
  geom_ribbon(data = filter(df, Model == "HER"), aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.2, fill = '#008080') +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey15", size = 0.1, width = 0.3) + 
  geom_line(aes(y = trillions, colour = Model,
                linetype = Model), size = 1.2) +  
  geom_point(aes(y = obs, shape = "Historical estimates from survey with 95% confidence intervals"),
             colour = "grey15", size = 2) +
  scale_shape_manual(values = 1) +
  theme(legend.position = c(0.35, 0.75)) +
  guides(shape  = guide_legend(order = 2),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 1)) +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  labs(x = "", y = NULL, shape = NULL) +
  ggtitle("Eggs spawned (trillions)") -> eggdep_plot
eggdep_plot
ggsave(paste0(fig_dir, "/compare_eggdep_CI.png"), plot = eggdep_plot, dpi = 300, height = 5, width = 10, units = "in")

# Posterior pred

df %>% 
  mutate(Model = factor(Model,
                        levels = c("HER", "LS"),
                        labels = c("HER with 95% credible & posterior predictive intervals", "LS"))) %>% 
  ggplot(aes(x = Year)) +
  geom_ribbon(data = filter(df, Model == "HER"), aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.3, fill = '#008080') +
  geom_ribbon(data = filter(df, Model == "HER"), aes(x = Year, ymin = p025, ymax = p975),
              alpha = 0.2, fill = '#008080') +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey15", size = 0.1, width = 0.3) + 
  geom_line(aes(y = trillions, colour = Model,
                linetype = Model), size = 1.2) +  
  geom_point(aes(y = obs, shape = "Historical estimates from survey with 95% confidence intervals"),
             colour = "grey15", size = 2) +
  scale_shape_manual(values = 1) +
  theme(legend.position = c(0.35, 0.75)) +
  guides(shape  = guide_legend(order = 2),
         color = guide_legend(order = 1),
         linetype = guide_legend(order = 1)) +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), breaks = seq(1980, YEAR-1, 5)) +
  labs(x = "", y = NULL, shape = NULL) +
  ggtitle("Eggs spawned (trillions)") -> eggdep_plot
eggdep_plot
ggsave(paste0(fig_dir, "/compare_eggdep_CI_PPI.png"), plot = eggdep_plot, dpi = 300, height = 5, width = 10, units = "in")

# Compare resids

data.frame(Year = D[["year"]],
           resids = D[["resd_egg_dep"]],
           Model = "HER") %>% 
  bind_rows(LS_byyear %>% 
              select(Year = year, resids = res_tot_egg) %>% 
              mutate(Model = "LS")) -> df

ggplot(df, aes(x = Year, y = resids, colour = Model, shape = Model)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point(size = 2) +
  facet_wrap(~ Model, ncol = 1) +
  labs(x = NULL, y = NULL) +
  ggtitle("Egg deposition residuals") +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), 
                     breaks = seq(1980, YEAR-1, 5)) +
  theme(legend.position = 'none') -> egg_resids
egg_resids
ggsave(paste0(fig_dir, "/egg_resids.png"), plot = egg_resids, dpi = 300, height = 5, width = 10, units = "in")



# Survival blocks ----
df <- data.frame(D[["year"]], 
                 D[["Mij"]])
colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
write.table(df,sep="\t", row.names=FALSE)

# Table summary
df %>%
  select(Year, M = `3`) %>%
  mutate(HER = exp(-M)) %>% 
  select(-M) %>% 
  left_join(LS_byage %>% 
              select(Year, LS = survival), by = "Year") %>% 
  group_by(LS, HER) %>% 
  mutate(Blocks = paste0(min(Year), "-", max(Year))) -> df_sum

df_sum %>% 
  distinct(Blocks, HER, LS) %>% 
  select(Blocks, HER, LS) -> survival_blks

survival_blks[c(2:3)] <- lapply(survival_blks[,c(2:3)], formatC, digits = 2, small.mark = ".", format = "f")

# Add in LS estimates
df %>%
  select(Year, M = `3`) %>% 
  mutate(survival = exp(-M),
         Model = "HER") %>% 
  bind_rows(LS_byage %>% 
              select(Year, survival) %>% 
              mutate(Model = "LS")) -> df_compare

# Figure summary
ggplot(df_compare, aes(x = Year, y = survival,  colour = Model, linetype = Model)) +
  geom_vline(xintercept = c(1998.5, 2014.5), colour = "lightgrey", linetype = 3, alpha = 0.4) +
  geom_line(size = 1) +
  # geom_point() +
  lims(y = c(0, 1)) +
  annotation_custom(tableGrob(survival_blks, rows = NULL, 
                              theme = ttheme_minimal(base_size = 14, base_colour = "black", 
                                                     parse = FALSE, 
                                                     padding = unit(c(4, 4), "mm"))), 
                    xmin = 1985, xmax = 2015, ymin = 0.1, ymax = 0.35) +
  ggtitle("Survival") +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), 
                     breaks = seq(1980, YEAR-1, 5)) +  
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(0.1, 0.8)) -> survival_plot
survival_plot
ggsave(paste0(fig_dir, "/survival.png"), plot = survival_plot, dpi = 300, height = 5, width = 10, units = "in")

# Plot model comparisons with 95% posterior interval on HER
df_sum %>% 
  ungroup() %>% 
  distinct() %>% 
  melt(id.vars = c("Year", "Blocks"), variable.name = "Model", value.name = "survival") %>%
  left_join(surv_sum %>% distinct(Blocks, q025, q975)) %>% 
  mutate(Model = factor(Model,
                        levels = c("HER", "LS"),
                        labels = c("HER with 95% credible interval", "LS"))) %>% 
  ggplot() +
  # annotation_custom(tableGrob(survival_blks, rows = NULL, 
  #                             theme = ttheme_minimal(base_size = 14, base_colour = "black", 
  #                                                    parse = FALSE, 
  #                                                    padding = unit(c(4, 4), "mm"))), 
  #                   xmin = 1985, xmax = 2015, ymin = 0.1, ymax = 0.35) +
  geom_ribbon(aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.2, fill = '#008080') +
  geom_line(aes(x = Year, y = survival, colour = Model, linetype = Model), size = 1) +
  # geom_vline(xintercept = c(1998.5, 2014.5), colour = "lightgrey", linetype = 3, alpha = 0.4) +
  lims(y = c(0, 1)) +
  ggtitle("Survival") +
  scale_x_continuous(limits = c(1980, YEAR-1), labels = seq(1980, YEAR-1, 5), 
                     breaks = seq(1980, YEAR-1, 5)) +  
  labs(x = NULL, y = NULL) +
  theme(legend.position = c(0.2, 0.8)) -> survival_plot2
survival_plot2
ggsave(paste0(fig_dir, "/survival_CI2.png"), plot = survival_plot2, dpi = 300, height = 5, width = 10, units = "in")

# Maturity/Selectivity ----

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
                      labels = c("3", "4", "5", "6", "7", "8+"))) %>% #,
         # To make sure selectivity is differentiable, it was scaled to have a
         # mean of 1 across all ages. This was done in log space by substracting
         # the mean from the vector of age-specific selectivities. See Tech Doc
         # p 11. Here we normalize it from 0 to 1.
         # proportion = ifelse(param == "Selectivity", (proportion - 0)/(max(proportion) - 0), 
         #                     proportion)) %>% 
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
  # geom_line(aes(linetype = `Time blocks`, colour = Model, group = `combos`)) +
  geom_line(aes(linetype = combos, colour = combos, group = combos), size = 1) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 3.5, colour = "grey", linetype = 2) +
  geom_line(aes(linetype = combos, colour = combos, group = combos), size = 1) +
  scale_color_manual(values = c('#008080', '#008080', '#DAA520', '#DAA520')) +
  scale_linetype_manual(values = c(1, 2, 1, 2)) +
  expand_limits(y = 0) +
  labs(x = "\nAge", y = "Proportion\n", linetype = NULL, colour = NULL) +
  ggtitle(paste0(par$param[1])) +
  theme(legend.position = c(0.7, 0.2),
        legend.key.size = unit(2,"line"),
        legend.spacing.y = unit(0, "cm"),
        plot.title = element_text(hjust = 0.5)) -> mat
mat
matsel %>% filter(param == "Selectivity") -> par

ggplot(par, aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = combos, colour = combos, group = combos), size = 1) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 3.5, colour = "grey", linetype = 2) +
  geom_line(aes(linetype = combos, colour = combos, group = combos), size = 1) +
  expand_limits(y = 0) +
  scale_color_manual(values = c('#008080', '#008080', '#DAA520')) +
  scale_linetype_manual(values = c(1, 2, 1)) +
  labs(x = "\nAge", y = NULL, linetype = NULL, color = NULL) +
  scale_y_continuous(breaks = seq(0, max(par$proportion), .25)) +
  ggtitle(paste0(par$param[1])) +
  theme(legend.position = c(0.7, 0.25),
        legend.key.size = unit(2,"line"),
        legend.spacing.y = unit(0, "cm"),
        plot.title = element_text(hjust = 0.5)) -> sel
sel
cowplot::plot_grid(mat, sel, align = "h", nrow = 1) -> matsel_plot

ggsave(paste0(fig_dir, "/mat_sel.png"), plot = matsel_plot, dpi = 300, height = 6, width = 10, units = "in")

# With uncertainty

mat_sum %>% 
  distinct(`Time blocks` = Blocks, Age, q025, q975) %>% 
  mutate(param = "Maturity") %>% 
  bind_rows(sel_sum %>% 
              distinct(`Time blocks` = Blocks, Age, q025, q975) %>% 
              mutate(param = "Selectivity")) %>% 
  mutate(age = as.numeric(as.character(Age))) %>% 
  select(-Age) %>% 
  right_join(matsel %>% 
              mutate(age = as.numeric(as.character(factor(Age, 
                                                          levels = c("3", "4", "5", "6", "7", "8+"),
                                                          labels = c("3", "4", "5", "6", "7", "8"))))),
            by = c("Time blocks", "param", "age")) -> matsel2

matsel2 %>% filter(param == "Maturity") -> par
ggplot() + 
  geom_line(data = par, 
            aes(x = age, y = proportion, linetype = combos, 
                colour = combos, group = combos), 
            size = 1) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 5.5, colour = "grey", linetype = 2) +
  geom_line(data = par, 
            aes(x = age, y = proportion, linetype = combos, 
                colour = combos, group = combos), 
            size = 1) +
  geom_ribbon(data = par %>% filter(Model == "HER") %>%
                distinct(age, `Time blocks`, q025, q975, param),
              aes(x = age, ymin = q025, ymax = q975, fill = `Time blocks`),
              alpha = 0.2) +
  scale_fill_manual(values = c('#008080', '#008080'), guide = FALSE) +
  scale_color_manual(values = c('#008080', '#008080', '#DAA520', '#DAA520')) +
  scale_linetype_manual(values = c(1, 2, 1, 2)) +
  expand_limits(y = 0) +
  labs(x = "\nAge", y = "Proportion\n", linetype = NULL, colour = NULL) +
  ggtitle(paste0(par$param[1])) +
  theme(legend.position = c(0.75, 0.2),
        legend.key.size = unit(2,"line"),
        legend.spacing.y = unit(0, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = c("3", "4", "5", "6", "7", "8+")) -> mat
mat

matsel2 %>% filter(param == "Selectivity") -> par

ggplot() + 
  geom_line(data = par, 
            aes(x = age, y = proportion, linetype = combos, 
                colour = combos, group = combos), 
            size = 1) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  geom_vline(xintercept = 5.5, colour = "grey", linetype = 2) +
  geom_line(data = par, 
            aes(x = age, y = proportion, linetype = combos, 
                colour = combos, group = combos), 
            size = 1) +
  geom_ribbon(data = par %>% filter(Model == "HER") %>%
                distinct(age, `Time blocks`, q025, q975, param),
              aes(x = age, ymin = q025, ymax = q975, fill = `Time blocks`),
              alpha = 0.2) +
  scale_fill_manual(values = c('#008080', '#008080'), guide = FALSE) +
  scale_color_manual(values = c('#008080', '#008080', '#DAA520', '#DAA520')) +
  scale_linetype_manual(values = c(1, 2, 1, 2)) +
  expand_limits(y = 0) +
  labs(x = "\nAge", y = NULL, linetype = NULL, colour = NULL) +
  ggtitle(paste0(par$param[1])) +
  theme(legend.position = c(0.7, 0.2),
        legend.key.size = unit(2,"line"),
        legend.spacing.y = unit(0, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = c("3", "4", "5", "6", "7", "8+")) -> sel
sel

cowplot::plot_grid(mat, sel, align = "h", nrow = 1) -> matsel_plot

ggsave(paste0(fig_dir, "/mat_sel_CI.png"), plot = matsel_plot, dpi = 300, height = 6, width = 10, units = "in")

# Age comps ----

df <- data.frame("Commercial fishery",
                 D[["data_cm_comp"]]) # change number of rows
colnames(df) <- c("Source", "Year", paste(D[['sage']]:D[['nage']]))
sp <- data.frame("Cast net survey",
                 D[["data_sp_comp"]]) # change number of rows
colnames(sp) <- c("Source", "Year", paste(D[['sage']]:D[['nage']]))
df <- bind_rows(df, sp); rm(sp)

df %>% 
  filter(Year >= 1980) %>% 
  melt(id.vars = c("Year", "Source"), variable.name = "Age", value.name = "obs") %>% 
  mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+"))) -> her_agecomps

# Model predictions and residuals
predcm <- data.frame("Commercial fishery",
                     D[["year"]],
                     D[["pred_cm_comp"]]) # change number of rows
colnames(predcm) <- c("Source", "Year", paste(D[['sage']]:D[['nage']]))

predsp <- data.frame("Cast net survey",
                     D[["year"]],
                     D[["pred_sp_comp"]]) # change number of rows
colnames(predsp) <- c("Source", "Year", paste(D[['sage']]:D[['nage']]))
df <- bind_rows(predcm, predsp); rm(predcm, predsp)

df %>% 
  melt(id.vars = c("Year", "Source"), variable.name = "Age", value.name = "pred") %>% 
  mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+"))) %>% 
  # join to observed age comps df
  left_join(her_agecomps) %>%
  group_by(Source) %>% 
  mutate(raw = obs - pred,
         # positive or negative
         `Model performance` = ifelse(raw >= 0, "Observed greater than estimated", "Observed less than estimated"),
         # sample size
         n = length(obs),
         # Raw
         resid = abs(raw),
         # Deviance residual http://data.princeton.edu/wws509/notes/c3s8.html,
         # deviance residual has the same direction as raw residual
         dev_resid = sqrt(2 * (obs * log(obs / pred ) + (n - obs) * log((n - obs)/(n - pred)))),
         direction = ifelse(raw >= 0, 1, -1),
         dev = direction * dev_resid,
         # Pearson's residuals
         pearson = (obs - pred)/ sqrt(var(pred))) -> her_agecomps

# 95% posterior predictive intervals
cm_comp_ppi <- cm_comp_ppi %>% 
  mutate(age = as.numeric(as.character(Age)),
         Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+")),
         Cohort = factor(ifelse(age < 8, Year - age, "plus"), ordered = TRUE))

# 95% posterior predictive intervals
sp_comp_ppi <- sp_comp_ppi %>% 
  mutate(age = as.numeric(as.character(Age)),
         Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+")),
         Cohort = factor(ifelse(age < 8, Year - age, "plus"), ordered = TRUE))

# Comp graph 1 ----
her_agecomps %>% # fishery
  filter(Source == "Cast net survey" &
           Year >= 2015) %>% 
  ggplot() + 
  geom_bar(aes(x = Age, y = obs), 
           colour = "white", fill = "lightgrey", stat = "identity", show.legend = FALSE,
           width = .9, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year, dir = "v", ncol = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Cast net survey")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) -> sp

#commercial catch age comps, purse seine
her_agecomps %>% # fishery
  filter(Source == "Commercial fishery"&
           Year >= 2015) %>% 
  ggplot() + 
  geom_bar(aes(x = Age, y = obs), 
           colour = "white", fill = "lightgrey", stat = "identity", show.legend = FALSE,
           width = .9, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year, dir = "v", ncol = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = NULL) + 
  ggtitle("Commercial fishery")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) -> cm

cowplot::plot_grid(sp, cm, align = "hv", nrow = 1) -> comp_plot
comp_plot
ggsave(paste0(fig_dir, "/comp_raw.png"), plot = comp_plot, dpi = 300, height = 7, width = 10, units = "in")

# Comp graph 2 ----

mycols <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "#a4a7b3" , "white")#BCBABE")
nages <- D[["nage"]] - D[["sage"]]
xtra <- nages * ((n_distinct(cm_comp_ppi$Cohort)-1)/nages - floor((n_distinct(cm_comp_ppi$Cohort)-1)/nages))
cohort_fill <- c(rep(mycols[1:5], (n_distinct(cm_comp_ppi$Cohort)-1)/nages), mycols[1:xtra], mycols[6])
cohort_cols <- replace(cohort_fill, cohort_fill == "white", "black") # Plus group always white bar with black border

df <- data.frame(spawnage_comp_obs = D[["data_sp_comp"]])
colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
df %>% 
  filter(Year >= 2015 ) %>%
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
              mutate(Model = "LS") ) %>% 
  filter(Year >= 2015) -> pred

her_agecomps %>% # fishery
  filter(Source == "Cast net survey" &
           Year >= 2015) %>% 
  left_join(sp_comp_ppi) %>%
  ggplot() + 
  geom_bar(aes(x = Age, y = obs, fill = Cohort), 
           colour = "white",stat = "identity", show.legend = FALSE,
           width = .9, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = cohort_fill) +
  facet_wrap(~ Year, dir = "v", ncol = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Cast net survey")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) -> sp

#commercial catch age comps, purse seine
df <- data.frame(comage_comp_obs = D[["data_cm_comp"]])
colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
df %>% 
  filter(Year >= 2015 ) %>%
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
              mutate(Model = "LS") ) %>% 
  filter(Year >= 2015) -> pred

her_agecomps %>% # fishery
  filter(Source == "Commercial fishery"&
           Year >= 2015) %>% 
  left_join(cm_comp_ppi) %>%
  ggplot() + 
  geom_bar(aes(x = Age, y = obs, fill = Cohort), 
           colour = "white",stat = "identity", show.legend = FALSE,
           width = .9, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = cohort_fill) +
  facet_wrap(~ Year, dir = "v", ncol = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = NULL) + 
  ggtitle("Commercial fishery")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) -> cm

cowplot::plot_grid(sp, cm, align = "hv", nrow = 1) -> comp_plot
comp_plot
ggsave(paste0(fig_dir, "/comp_cohort.png"), plot = comp_plot, dpi = 300, height = 7, width = 10, units = "in")

# Comp graph 3 ----

df <- data.frame(spawnage_comp_obs = D[["data_sp_comp"]])
colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
df %>% 
  filter(Year >= 2015 ) %>%
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
              mutate(Model = "LS") ) %>% 
  filter(Year >= 2015) -> pred

her_agecomps %>% # fishery
  filter(Source == "Cast net survey" &
           Year >= 2015) %>% 
  left_join(sp_comp_ppi) %>%
  ggplot() + 
  geom_bar(aes(x = Age, y = obs, fill = Cohort), 
           colour = "white",stat = "identity", show.legend = FALSE,
           width = .9, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = cohort_fill) +
  geom_line(data = pred, aes(x = Age, y = proportion, colour = Model, 
                             linetype = Model, group = Model), size = 1.5) +
  facet_wrap(~ Year, dir = "v", ncol = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Cast net survey")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) -> sp

#commercial catch age comps, purse seine
df <- data.frame(comage_comp_obs = D[["data_cm_comp"]])
colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
df %>% 
  filter(Year >= 2015 ) %>%
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
              mutate(Model = "LS") ) %>% 
  filter(Year >= 2015) -> pred

her_agecomps %>% # fishery
  filter(Source == "Commercial fishery"&
           Year >= 2015) %>% 
  left_join(cm_comp_ppi) %>%
  ggplot() + 
  geom_bar(aes(x = Age, y = obs, fill = Cohort), 
           colour = "white",stat = "identity", show.legend = FALSE,
           width = .9, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = cohort_fill) +
  geom_line(data = pred, aes(x = Age, y = proportion, colour = Model, 
                             linetype = Model, group = Model), size = 1.5) +
  facet_wrap(~ Year, dir = "v", ncol = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = NULL) + 
  ggtitle("Commercial fishery")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) -> cm

cowplot::plot_grid(sp, cm, align = "hv", nrow = 1) -> comp_plot
comp_plot
ggsave(paste0(fig_dir, "/comp_pred.png"), plot = comp_plot, dpi = 300, height = 7, width = 10, units = "in")

# Comp graph 4 ----

df <- data.frame(spawnage_comp_obs = D[["data_sp_comp"]])
colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
df %>% 
  filter(Year >= 2015 ) %>%
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
              mutate(Model = "LS") ) %>% 
  filter(Year >= 2015) -> pred

her_agecomps %>% # fishery
  filter(Source == "Cast net survey" &
           Year >= 2015) %>% 
  left_join(sp_comp_ppi) %>%
  ggplot() + 
  geom_bar(aes(x = Age, y = obs, fill = Cohort), 
           colour = "white",stat = "identity", show.legend = FALSE,
           width = .9, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = cohort_fill) +
  geom_errorbar(aes(x = Age, ymin = q025, ymax = q975),
                colour = "#008080", width = .1, size = 1, position = position_dodge(width = 0.5)) +
  geom_line(data = pred, aes(x = Age, y = proportion, colour = Model, 
                             linetype = Model, group = Model), size = 1.5) +
  facet_wrap(~ Year, dir = "v", ncol = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Cast net survey")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) -> sp

#commercial catch age comps, purse seine
df <- data.frame(comage_comp_obs = D[["data_cm_comp"]])
colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
df %>% 
  filter(Year >= 2015 ) %>%
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
              mutate(Model = "LS") ) %>% 
  filter(Year >= 2015) -> pred

her_agecomps %>% # fishery
  filter(Source == "Commercial fishery"&
           Year >= 2015) %>% 
  left_join(cm_comp_ppi) %>%
  ggplot() + 
  geom_bar(aes(x = Age, y = obs, fill = Cohort), 
           colour = "white",stat = "identity", show.legend = FALSE,
           width = .9, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = cohort_fill) +
  geom_errorbar(aes(x = Age, ymin = q025, ymax = q975),
                colour = "#008080", width = .1, size = 1, position = position_dodge(width = 0.5)) +
  geom_line(data = pred, aes(x = Age, y = proportion, colour = Model, 
                             linetype = Model, group = Model), size = 1.5) +
  facet_wrap(~ Year, dir = "v", ncol = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = NULL) + 
  ggtitle("Commercial fishery")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) -> cm

cowplot::plot_grid(sp, cm, align = "hv", nrow = 1) -> comp_plot
ggsave(paste0(fig_dir, "/comp_CI.png"), plot = comp_plot, dpi = 300, height = 7, width = 10, units = "in")

# Weight-at-age ----

df <- as.data.frame(D[["data_cm_waa"]])
colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
df <-	df %>% 
  gather("Age", "weight", -Year) %>% 
  mutate(age = as.integer(Age),
         Cohort = as.factor(Year - age))

# 1) by cohort

# Create a better colour palette (lots of options, see helper.r library ggthemes)
pal <- ggthemes::canva_pal("Warm and cool")(4) 

# Axis ticks for plot (see helper.r tickr() fxn for details)
axis <- tickr(df, Year, 5)
ggplot(df, aes(Year, weight, colour = Cohort)) +
  geom_line(size = 1) +#alpha = 0.90) +
  geom_point(aes(fill = Cohort), show.legend = FALSE, size = 1) +
  labs(x = "Year", y = "Weight-at-age (grams)\n", colour = "Cohort") +
  guides(colour = guide_legend(ncol = 9)) +
  scale_colour_manual(values = colorRampPalette(pal)(n_distinct(df$Cohort))) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) -> waa_cohort_plot

ggsave(paste0(fig_dir, "/waa_cohort_plot.png"), plot = waa_cohort_plot, dpi = 300, height = 4, width = 6, units = "in")

# 2) by age
df %>% 
  group_by(Age) %>% 
  summarize(mean_weight = mean(weight)) -> means

ggplot(df, aes(Year, weight, group = Age, colour = Age)) + 
  geom_line() + 
  geom_point() +
  geom_hline(data = means, aes(colour = Age, yintercept = mean_weight), alpha = 0.4, linetype = 2) + 
  labs(x = "Year", y = "Weight-at-age (grams)\n", colour = "Age") +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = colorRampPalette(pal)(n_distinct(df$Age))) +
  guides(colour = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) -> waa_plot

ggsave(paste0(fig_dir, "/waa_plot.png"), plot = waa_plot, dpi = 300, height = 4, width = 6, units = "in") 

# Derive SR params ----

# I wrote this code to derive phie and thus so and beta parameters from the stock
# recruitment relationship because I was unable to pull beta from D via
# D[["beta"]] (rounding issue?)

# df <- as.data.frame(D[["data_sp_waa"]])
# colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
# df %>% 
#   gather("Age", "weight", -Year) %>% 
#   filter(Year >= D[["mod_syr"]]) %>% 
#   group_by(Age) %>% 
#   summarise(waa = mean(weight)) -> waa_bar
# 
# df <- as.data.frame(D[["mat"]])
# colnames(df) <- c(paste(D[['sage']] : D[['nage']]))
# df %>% 
#   mutate(Year = D[["year"]]) %>% 
#   gather("Age", "maturity", -Year) %>% 
#   filter(Year >= D[["mod_syr"]]) %>% 
#   group_by(Age) %>% 
#   summarise(mat = mean(maturity)) -> mat_bar
# 
# df <- as.data.frame(D[["Mij"]])
# colnames(df) <- c(paste(D[['sage']] : D[['nage']]))
# df %>% 
#   mutate(Year = D[["year"]]) %>% 
#   gather("Age", "mort", -Year) %>% 
#   filter(Year >= D[["mod_syr"]]) %>% 
#   #group_by(Age) %>% 
#   summarise(mort = mean(mort)) %>% 
#   pull(mort) -> mort_bar
# 
# lx <- c(1:6)
# 
# lx[1] <- 1
# 
# for(i in 2:6){
#   lx[i] <- lx[i-1] * exp(-mort_bar)
# }
# # plus group
# lx[6] = lx[6] / (1 - exp(-mort_bar))
# lx  
# 
# waa_bar %>% 
#   left_join(mat_bar) %>% 
#   mutate(lx = lx,
#          fec = waa * mat,
#          phi = lx * fec) %>% 
#   summarize(phie = sum(phi)) %>% 
#   pull(phie) -> phie
# 
# # From .par file 
# log_ro <- D[["theta"]][4] # theta(4)
# log_reck <- D[["theta"]][5]  # theta(5)
# 
# ro <- exp(log_ro) # equilibrium recruitment
# reck <-  exp(log_reck) + 1.0 # recruitment compensation ratio
# 
# # stock-recruit parameters
# so <- reck/phie;
# beta <- log(reck) / (ro * phie) # beta = (reck - 1.0)/(ro *phiE)
# 
# so; beta
# D[["so"]];D[["beta"]]

# her_test.r code ----
# 
# # HER - ADFG's herring model. Original ADMB code written by SJD Martell. Helper
# # files and documention contributed by M. Rudd.
# 
# source("R/tools.R") 
# library(R2admb)
# 
# 
# # clean up directory - remove unnecessary files to run
# # add ".pin" when wanting to save .pin file, or any other file to save
# setwd("HER/")
# need_files <- c(".tpl", ".ctl", ".dat", ".R")
# files_present <- list.files()
# keep_index <- unlist(sapply(1:length(need_files), function(x) which(grepl(need_files[x], files_present))))
# rm_files <- files_present[-keep_index]
# remove <- sapply(1:length(rm_files), function(x) unlink(rm_files[x], TRUE))
# # setwd("..") # backs out one folder
# 
# ## compile model
# compile_admb("her", verbose=TRUE)
# ## run MAP
# run_admb("her")
# 
# years <- 1971:2015
# 
# #Maturity
# readMat("mat", file="her.rep", nrow = length(years))
# readVec("mat_params[1]", file="her.rep")
# readVec("mat_params[2]", file="her.rep")
# 
# # Natural mortality
# natmat <- readMat("Mij", file = "her.rep", nrow = length(years))[,1]
# plot(natmat ~ years, ylim=c(0, max(natmat)*1.1), type="l", lwd=2)
# 
# ## read report from initial MAP run - maximum a posteriori estimation (i.e.
# ## maximum likelihood using priors!)
# ssb <- readVec("ssb", file="her.rep")
# ## put in thousands
# ssb <- ssb/1000
# 
# plot(ssb ~ years, ylim=c(0, max(ssb)*1.1), type="l", lwd=2)
# 
# 
# ## run simulation with seed 123
# run_admb("her", extra.args="-sim 123")
# 
# ## run MCMC
# run_admb("her", extra.args="-mcmc 10000 -mcsave 10")
# run_admb("her", extra.args="-mceval")
# 
# ## posterior distributions
# ssb_ps <- read.table("ssb.ps")
# natural_ps <- read.table("natural.ps", header=TRUE)
# 
# par(mfrow=c(2,1))
# plot(natural_ps[,1])
# abline(h=median(natural_ps[,1]), col="red")
# plot(natural_ps[,2])
# abline(h=median(natural_ps[,2]), col="red")
# 
# ## from tech doc:
# ## 1) first fit to sitka data
# ## 2) save .par file as her.pin
