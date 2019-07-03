# Run HER and compare results to LS
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-06-03

# HER - ADFG's new herring model. Original ADMB code written by SJD Martell. Helper
# files and documentation contributed by M. Rudd.

# This code creates figures for the least-squares (LS) model and HER and a set
# of comparison figures.

# Set up ----

set.seed(907)

library(ggdistribute) # plot posteriors
library(coda) # mcmc diagnostics
library(BayesianTools) # more mcmc diagnostics

# Forecast year
YEAR <- 2019

# User input: do you want to run the full MCMC?
run_mcmc <- TRUE

# Create directory for HER/LS comparison figs, HER figs, and LS figs
main_dir <- getwd()
# Directory for all comparison figures
fig_dir <- file.path(main_dir, paste0(YEAR, "_forecast/results/compare_HER_LS"))
dir.create(fig_dir, showWarnings = FALSE)
# Path for HER-specific figures
HERfig_dir <- file.path(fig_dir, "HER")
dir.create(HERfig_dir, showWarnings = FALSE)
# Subdirectory of HER MCMC diagnostic figures
HERmcmc_dir <- file.path(HERfig_dir, "MCMC_diagnostics")
dir.create(HERmcmc_dir, showWarnings = FALSE)
# Path for LS-specific figures
LSfig_dir <- file.path(fig_dir, "LS")
dir.create(LSfig_dir, showWarnings = FALSE)

# source(paste0(YEAR, "_forecast/r/tools.r"))
source(paste0(YEAR, "_forecast/r/helper.r"))

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
setwd(file.path(main_dir, paste0(YEAR, "_forecast/admb/HER_bestLS")))

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

# Run MCMC ----

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

# NOTE: See helper.R for source code for all user-defined functions.

# Posterior intervals (aka credible intervals, show variability around
# the expected value)
egg_sum <- ps_byyear(fn = "pred_egg_dep")
spb_sum <- ps_byyear(fn = "sp_B", unit_conv = 0.90718)
matb_sum <- ps_byyear(fn = "mat_B", unit_conv = 0.90718)
catch_sum <- ps_byyear(fn = "pred_catch") # warning msg ok!
ricker_sum <- ps_byyear(fn = "ricker_rec", syr = D[["mod_syr"]] + D[["sage"]],
                         lyr = D[["mod_nyr"]] + 1)
age3_sum <- ps_byyear(fn = "age3_rec", syr = D[["mod_syr"]],
                         lyr = D[["mod_nyr"]] + 1)
fmort_sum <- ps_byyear(fn = "fmort", syr = D[["mod_syr"]], # won't be anything if conditioned on catch
                         lyr = D[["mod_nyr"]])

# Survival posterior distribution
surv_sum <- ps_byage(fn = "survival", syr = D[["mod_syr"]],
                    lyr = D[["mod_nyr"]], n = niter / thin + 1)

# Maturity posterior distribution
mat_sum <- ps_byage(fn = "maturity", syr = D[["mod_syr"]], 
                    lyr = D[["mod_nyr"]], n = niter / thin + 1)

# Selectivity posterior distribution
sel_sum <- ps_byage(fn = "selectivity", syr = D[["mod_syr"]], 
                    lyr = D[["mod_nyr"]] + 1, n = niter / thin + 1)  
sel_sum <- sel_sum[min != YEAR]

# Posterior intervals (credible intervals) for age compositions
sp_comp_sum <- ps_comps(fn = "pred_sp_comp", syr = D[["mod_syr"]], 
                        lyr = D[["mod_nyr"]], n = niter / thin + 1) 

cm_comp_sum <- ps_comps(fn = "pred_cm_comp", syr = D[["mod_syr"]], 
                    lyr = D[["mod_nyr"]], n = niter / thin + 1)  

# Posterior predictive intervals 

# Variance estimates for each posterior sample from the multivariate logistic
# distribution used for age comps
sp_tau2 <- ps_tau(fn = "pp_sp_tau2")
cm_tau2 <- ps_tau(fn = "pp_cm_tau2")

# Get posterior predictive intervals for multivariate logistic distribution (age
# compositions). Uses the posterior sample estimates and variance.
cm_comp_ppi <- ppi_rmvlogistic(ps_sum = cm_comp_sum,
                               tau2 = cm_tau2)

sp_comp_ppi <- ppi_rmvlogistic(ps_sum = sp_comp_sum,
                               tau2 = sp_tau2)

# MCMC diagnostic plots ----

# Parameter draws
pars <- D[["post.samp"]]
burn <- round(burn_in * (niter / thin + 1), 0) # index for end of burn in
pars <- pars[burn:(niter/thin),] # Remove burn-in

ages <- c(D[["sage"]]:D[["nage"]])
yrs <- c(D[["mod_syr"]]:D[["mod_nyr"]])

# Labels for maturity blocks
out <- list()
mat_sum <- mat_sum %>% filter(Blocks != "1980-2018")
for(i in 1:length(unique(mat_sum$Blocks))) {
  out[[i]] <- c(paste0("mat_a50_", i), paste0("mat_a95_", i))
}

mat_labels <- unlist(out)

# Labels for selectivity blocks
out <- list()
for(i in 1:length(unique(sel_sum$Blocks))) {
  out[[i]] <- c(paste0("log_s50_", i), paste0("log_selk_", i))
}
sel_labels <- unlist(out)

colnames(pars) <- c("log_Mbar", "log_rinit", "log_rbar", "log_ro", "log_reck", 
                    paste0("log_rinit_devs_", ages), paste0("log_rbar_devs_", yrs),
                    mat_labels, paste0("log_Mdevs_", 1:length(unique(surv_sum$Blocks))),
                    sel_labels)

pars <- as.data.frame(pars) %>% 
  select(- contains("rinit_devs"), - contains("rbar_devs")) 

pars <- coda::as.mcmc(pars)

# Autocorrelation figures (ACF): There shouldn't be any correlation in lag >= 1.
# ACFs look great except for log_rinit, although the correlation declines
# steeply, so I don't think this is too much of a concern.
acfplot(pars)
png(paste0(HERmcmc_dir, "/acf.png"))
acfplot(pars)
dev.off()

# Trace plots (aka caterpillar plots) with posterior marginal densities:
# Trace plots shouldn't have any trends, should be well-mixed.
png(paste0(HERmcmc_dir, "/param_trace.png"),
    width = 11, height = 11, units = "in", res = 300)
par(mfrow = c(7, 4))
plot(pars, auto.layout = FALSE, ask = FALSE)
dev.off()

# Marginal densities (diagonal), pairwise densities (lower panes) and
# correlation coefficien (upper panels) - when there is strong correlation, the
# marginal density (uncertainty) in the parameter estimate is affected.
png(paste0(HERmcmc_dir, "/param_correlation.png"),
    width = 11, height = 11, units = "in", res = 300)
BayesianTools::correlationPlot(data.frame(pars)) 
dev.off()

# geweke.diag returns Z scores for the equality of (by default) the first 10%
# and the last 50% of the chain.  pnorm(abs(v),lower.tail=FALSE)*2 computes a
# two-tailed Z-test, of which the results should be > 0.05 if chain is converged)
pnorm(geweke.diag(pars)$z,lower.tail=FALSE)*2
png(paste0(HERmcmc_dir, "/geweke_diagnostic.png"),
    width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(5, 3))
geweke.plot(pars, auto.layout = FALSE) # should stay within the 2 sd bounds
dev.off()
# More on Geweke:
# https://www2.math.su.se/matstat/reports/master/2011/rep2/report.pdf or
# http://pymc-devs.github.io/pymc/modelchecking.html?highlight=geweke

# Compute the effective sample size, corrected for autocorrelation. This
# value should be > 200 for reasonal estimation of credible intervals
coda::effectiveSize(pars)

# Survival is programmed such that it can be age-dependent but is currently
# constant across all ages.
surv_sum %>% 
  filter(Age == "3") %>% 
  mutate(Age = "All ages") -> surv_sum
mcmc_plot(df = surv_sum, type = "ps_byage", name = "Survival")
mcmc_plot(df = mat_sum, type = "ps_byage", name = "Maturity")
mcmc_plot(df = sel_sum, type = "ps_byage", name = "Selectivity")

# Derived posterior trajectories
mcmc_plot(df = spb_sum, type = "ps_byyear", "Spawning biomass", height = 10, width = 12)
mcmc_plot(df = matb_sum, type = "ps_byyear", "Mature biomass", height = 10, width = 12)
mcmc_plot(df = egg_sum, type = "ps_byyear", "Egg deposition", height = 10, width = 12)
mcmc_plot(df = catch_sum, type = "ps_byyear", "Predicted catch", height = 10, width = 12)
mcmc_plot(df = ricker_sum, type = "ps_byyear", "Ricker-estimated age-3 recruitment", height = 10, width = 12)
mcmc_plot(df = age3_sum, type = "ps_byyear", "ASA-estimated age-3 recruitment", height = 10, width = 12)
# mcmc_plot(df = fmort_sum, type = "ps_byyear", "Fishing mortality") # If conditioned on catch
mcmc_plot(df = cm_comp_sum, type = "ps_comps", "Commercial fishery age composition", height = 20, width = 14)
mcmc_plot(df = sp_comp_sum, type = "ps_comps", "Cast net survey age composition", height = 20, width = 14)

# Posterior distribution of forecasted mature biomass 
df <- data.frame(read.table("fore_matb.ps"))
colnames(df) <- "quantity"
df %>% 
  mutate(Year = D[["mod_nyr"]] + 1,
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

fore_matb_sum %>% 
  ggplot(aes(x = iter, y = quantity)) +
  geom_hline(aes(yintercept = mean), linetype = 2, colour = "grey") +
  geom_line() -> caterpillar_forematb

ggsave(paste0(HERmcmc_dir, "/caterpillar_forematb.png"), plot = caterpillar_forematb, dpi = 300, height = 4, width = 6, units = "in")

fore_matb_sum %>% 
  ggplot(aes(x = quantity)) +
  geom_posterior() +
  labs(x = "\nForecast mature biomass posterior distribution") + 
  scale_x_continuous(labels = scales::comma) -> post_forematb

ggsave(paste0(HERmcmc_dir, "/posterior_forematb.png"), plot = post_forematb, dpi = 300, height = 4, width = 6, units = "in")

# Diagnostics ----
P <- read_fit("her")
P[["nopar"]]
P[["nlogl"]]
P[["logDetHess"]]
P[["maxgrad"]]

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
  mutate(Model = "HER") %>% 
  bind_rows(
    LS_byyear %>% 
      select(Year = year, mean = tot_mat_B_tons) %>% 
      mutate(type = "Estimate") %>% 
      bind_rows(
        LS_forec %>% 
          summarize(mean = sum(for_mat_baa_tons)) %>% 
          mutate(Year = max(LS_byyear$year) + 1,
                 type = "Forecast")) %>% 
      mutate(Model = "LS")) %>% 
  mutate(label = ifelse(type == "Forecast", 
                        prettyNum(mean, big.mark = ",", digits = 1), NA)) -> df
tickryr <- data.frame(Year = 1980:2020)
axisf <- tickr(tickryr, Year, 5)
ggplot(df, aes(x = Year)) +
  geom_ribbon(data = filter(df, !is.na(q025)), aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.4, fill = "grey80") +
  geom_line(aes(y = mean, colour = Model, linetype = Model, group = Model), size = 1) +
  geom_point(data = LS_byyear, aes(x = year, y = surv_est_matbio, shape = "Historical estimates from survey")) +
  geom_point(data = df %>% filter(type == "Forecast"), aes(x = Year, y = mean, colour = Model), shape = 8) +
  geom_text_repel(data = df, aes(x = Year, y = mean, label = label, colour = Model),
                  size = 3, nudge_x = 1.5, show.legend = FALSE) +
  scale_shape_manual(values = 1) +
  scale_colour_grey() +
  scale_x_continuous(limits = c(min(tickryr$Year), max(tickryr$Year)),
                     breaks = axisf$breaks, labels = axisf$labels) +
  scale_y_continuous(limits = c(0, max(LS_byyear$surv_est_matbio)), 
                     labels = scales::comma) +
  labs(x = NULL, y = "Mature biomass (tons)\n", shape = NULL, text = NULL) +
  theme(legend.position = c(0.25, 0.8)) -> post_matb

ggsave(paste0(fig_dir, "/compare_matbiomass.png"), plot = post_matb, dpi = 300, height = 4, width = 6, units = "in")

# For HER with posterior distributions: 5%, 50%, and 95% credibility intervals
# (or posterior predictive interval)

# Combine estimates and forecast
matb_sum %>% 
  mutate(Year = as.numeric(as.character(Year)),
         type = "Estimate") %>% 
  bind_rows(fore_matb_sum %>% 
              mutate(type = "Forecast")) %>% 
  distinct(Year, mean, q025, q975, q250, q750, q475, q525, type) %>% 
  mutate(label = ifelse(type == "Forecast", 
                        prettyNum(mean, big.mark = ",", digits = 1), NA)) -> df

tickryr <- data.frame(Year = 1980:YEAR+3)
axisf <- tickr(tickryr, Year, 5)

ggplot(df, aes(x = Year)) +
  geom_ribbon(aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.6, fill = "grey70") +
  geom_ribbon(aes(x = Year, ymin = q250, ymax = q750),
              alpha = 0.6, fill = "grey40") +
  geom_ribbon(aes(x = Year, ymin = q475, ymax = q525),
              alpha = 0.6, fill = "black") +
  geom_point(data = LS_byyear, aes(x = year, y = surv_est_matbio, shape = "Historical estimates from survey")) +
  geom_point(data = df %>% filter(type == "Forecast"), aes(x = Year, y = mean), shape = 8) +
  geom_text_repel(data = df, aes(x = Year, y = mean, label = label),
                  size = 3, nudge_x = 1.5, show.legend = FALSE) +
  scale_shape_manual(values = 1) + 
  scale_colour_grey() +
  scale_x_continuous(limits = c(min(tickryr$Year), max(tickryr$Year)),
                     breaks = axisf$breaks, labels = axisf$labels) +
  scale_y_continuous(limits = c(0, max(LS_byyear$surv_est_matbio)), labels = scales::comma) +
  labs(x = NULL, y = "Mature biomass (tons)\n", shape = NULL, text = NULL) +
  theme(legend.position = c(0.25, 0.8)) -> post_matb

ggsave(paste0(HERfig_dir, "/matbiomass_CI.png"), plot = post_matb, dpi = 300, height = 4, width = 6, units = "in")

# Spawning biomass ----

# post-fishery
spb_sum %>% 
  mutate(Year = as.numeric(as.character(Year)),
         type = "Estimate") %>% 
  distinct(Year, mean, q025, q975, type) %>% 
  mutate(Model = "HER") %>% 
  bind_rows(
    LS_byyear %>% 
      select(Year = year, mean = tot_sp_B_tons) %>% 
      mutate(type = "Estimate",
             Model = "LS")) -> df

srv_index <- data.frame(Year = LS_byyear$year,
                        srv_spB = LS_byyear$surv_est_spB,
                        Model = "Historical survey index")

axisx <- tickr(df, Year, 5)
ggplot() +
  geom_ribbon(data = filter(df, !is.na(q025)), aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.4, fill = "grey80") +
  geom_line(data = df, aes(x = Year, y = mean, colour = Model, linetype = Model), size = 1) +
  geom_point(data = srv_index, aes(x = Year, y = srv_spB, shape = "Historical estimates from survey")) +
  scale_shape_manual(values = 1) +
  scale_colour_grey() +
  theme(legend.position = c(0.25, 0.7)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  scale_y_continuous(limits = c(0, max(LS_byyear$surv_est_spB)), labels = scales::comma) +
  labs(x = "", y = "Spawning biomass (tons)\n", shape = "Data") -> spbiomass_plot

ggsave(paste0(fig_dir, "/compare_spbiomass_plot.png"), plot = spbiomass_plot, dpi = 300, height = 4, width = 6, units = "in")

# For HER with credibility intervals:
spb_sum %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  distinct(Year, mean, q025, q975, q250, q750, q475, q525) -> df

axisx <- tickr(df, Year, 5)
ggplot(data = df, aes(x = Year)) +  
  geom_ribbon(aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.6, fill = "grey70") +
  geom_ribbon(aes(x = Year, ymin = q250, ymax = q750),
              alpha = 0.6, fill = "grey40") +
  geom_ribbon(aes(x = Year, ymin = q475, ymax = q525),
              alpha = 0.6, fill = "black") +
  geom_point(data = LS_byyear, aes(x = year, y = surv_est_spB, shape = "Historical estimates from survey")) +
  scale_colour_grey() +
  scale_shape_manual(values = 1) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  scale_y_continuous(limits = c(0, max(LS_byyear$surv_est_spB)), labels = scales::comma) +
  labs(x = NULL, y = "Spawning biomass (tons)\n", shape = NULL) +
  theme(legend.position = c(0.25, 0.8)) -> post_ssb

ggsave(paste0(HERfig_dir, "/spbiomass_CI.png"), plot = post_ssb, dpi = 300, height = 4, width = 6, units = "in")

# LS past model comparison ----

# bind past estimates/forecast of mature biomass with current
LS_yearminus %>% 
  bind_rows(LS_byyear %>% 
              select(Year = year, matbio_tons = tot_mat_B_tons) %>% 
              mutate(type = "estimate",
                     year_minus = 0),
            # current year forecast
            LS_forec %>% 
              summarize(matbio_tons = sum(for_mat_baa_tons)) %>% 
              mutate(Year = max(LS_byyear$year) + 1,
                     type = "forecast",
                     year_minus = 0)) %>% 
  group_by(year_minus) %>% 
  mutate(`Forecast model` = as.character(max(Year))) %>% 
  ungroup() %>% 
  mutate(`Forecast model` = fct_rev(factor(`Forecast model`)),
         label = ifelse(type == "forecast", 
                        prettyNum(matbio_tons, big.mark = ",", digits = 1), NA)) -> df
tickryr <- data.frame(Year = 1980:YEAR+3)
axisf <- tickr(tickryr, Year, 5)
ggplot() +
  geom_line(data = df, aes(x = Year, y = matbio_tons, 
                           colour = `Forecast model`, linetype = `Forecast model`)) +
  geom_point(data = df %>% filter(type == "forecast"), 
             aes(x = Year, y = matbio_tons, colour = `Forecast model`)) +
  geom_point(data = LS_byyear, 
             aes(x = year, y = surv_est_matbio,   shape = "Historical estimates from survey")) +
  geom_text_repel(data = df, aes(x = Year, y = matbio_tons, label = label, colour = `Forecast model`),
                  size = 2) +
  scale_shape_manual(values = 1) +
  scale_colour_grey() +
  theme(legend.position = c(0.25, 0.8),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(limits = c(min(tickryr$Year), max(tickryr$Year)),
                     breaks = axisf$breaks, labels = axisf$labels) +
  scale_y_continuous(limits = c(0, max(LS_byyear$surv_est_matbio)),
                     labels = scales::comma) +
  labs(x = "", y = "Mature biomass (tons)\n", shape = NULL) -> past_matbio

ggsave(paste0(LSfig_dir, "/compare_past_matbio.png"), plot = past_matbio, dpi = 300, height = 4, width = 6, units = "in")

# Biom/Abd bargraphs ----

# For LS:

LS_byyear %>%
  select(year, tcb, tot_sp_B_tons, N, tot_sp_N, tot_mat_N) %>% 
  mutate(tcn = tot_mat_N - tot_sp_N) %>% 
  gather("var", "value", -c(year)) -> df

axisb <- tickr(df, year, 5)

# Mature biomass bargraph
ggplot(df %>% 
         filter(var %in% c("tot_sp_B_tons", "tcb")) %>% 
         mutate(var = factor(var, ordered = TRUE, 
                             levels = c("tot_sp_B_tons", "tcb"),
                             labels = c("Spawning biomass", "Catch"))),
       aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  scale_fill_manual(values = c("white", "grey70")) +
  theme(legend.position = c(0.15, 0.85),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma, 
                     limits = c(0, max(LS_byyear$tot_sp_B_tons) * 1.2)) +
  labs(x = NULL, y = "Mature biomass (tons)", linetype = NULL, fill = NULL) -> matb_bar

# Mature abundance bargraph
ggplot(df %>% 
         filter(var %in% c("tot_sp_N", "tcn")) %>% 
         mutate(var = factor(var, ordered = TRUE, 
                             levels = c("tot_sp_N", "tcn"),
                             labels = c("Spawning abundance", "Catch"))),
       aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  scale_fill_manual(values = c("white", "grey70")) +
  theme(legend.position = c(0.15, 0.8),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma, 
                     limits = c(0, max(LS_byyear$tot_sp_N) * 1.65)) +
  labs(x = NULL, y = "Mature abundance (millions)", linetype = NULL, fill = NULL) -> matn_bar

# Total abundance bargraph
ggplot(df %>% 
         filter(var %in% c("N", "tcn")) %>% 
         mutate(var = factor(var, ordered = TRUE, 
                             levels = c("N", "tcn"),
                             labels = c("Mature + immature abundance", "Catch"))),
       aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  scale_fill_manual(values = c("white", "grey70")) +
  theme(legend.position = c(0.2, 0.85),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma, 
                     limits = c(0, max(LS_byyear$N) * 1.5)) +
  labs(x = NULL, y = "Total abundance (millions)", linetype = NULL, fill = NULL) -> totn_bar

cowplot::plot_grid(matb_bar, matn_bar, totn_bar, align = "hv", nrow = 3) -> bars

ggsave(paste0(LSfig_dir, "/biom_abd_barplots.png"), plot = bars, dpi = 300, height = 8, width = 7, units = "in")

# For HER:

tot_yrs <- D[["dat_nyr"]] - D[["dat_syr"]] + 1 # includes data
mod_yrs <- D[["mod_nyr"]] - D[["mod_syr"]] + 1 # just model years

# Spawning abundance: Convert sp_Nij to sp_Ni in millions
df <- D[["sp_Nij"]][1:mod_yrs,] # exclude the forecast year
df <- as.data.frame(df) 
colnames(df) <- c(paste(D[['sage']]:D[['nage']]))

sp_N <- df %>% mutate(year = (D[["mod_syr"]]):(D[["mod_nyr"]])) %>% 
  gather("Age", "sp_N", -year) %>% 
  group_by(year) %>% 
  summarize(sp_N = sum(sp_N))

# Catch in numbers: Convert Cij to Ci in millions
df <- D[["Cij"]][1:mod_yrs,] # exclude the forecast year
df <- as.data.frame(df) 
colnames(df) <- c(paste(D[['sage']]:D[['nage']]))

C <- df %>% mutate(year = (D[["mod_syr"]]):(D[["mod_nyr"]])) %>% 
  gather("Age", "C", -year) %>% 
  group_by(year) %>% 
  summarize(C = sum(C))

# Total abundance: Convert Nij to Ni in millions
df <- D[["Nij"]][1:mod_yrs,] # exclude the forecast year
df <- as.data.frame(df) 
colnames(df) <- c(paste(D[['sage']]:D[['nage']]))

N <- df %>% mutate(year = (D[["mod_syr"]]):(D[["mod_nyr"]])) %>% 
  gather("Age", "N", -year) %>% 
  group_by(year) %>% 
  summarize(N = sum(N))

# Spawning biomass and catch biomass
sp_B <- data.frame(year = (D[["mod_syr"]]):(D[["mod_nyr"]]),
                 sp_B = D[["sp_B"]] / 0.90718,
                 catch = D[["data_catch"]][10:tot_yrs, 2]  / 0.90718)

# Join and gather indices
sp_B %>% left_join(sp_N) %>% left_join(C) %>% left_join(N) %>% 
  gather("var", "value", -c(year)) -> df

# Mature biomass bargraph
ggplot(df %>% 
         filter(var %in% c("sp_B", "catch")) %>% 
         mutate(var = factor(var, ordered = TRUE, 
                             levels = c("sp_B", "catch"),
                             labels = c("Spawning biomass", "Catch"))),
       aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  scale_fill_manual(values = c("white", "grey70")) +
  theme(legend.position = c(0.15, 0.85),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma, 
                     limits = c(0, max(sp_B$sp_B) * 1.2)) +
  labs(x = NULL, y = "Mature biomass (tons)", linetype = NULL, fill = NULL) -> matb_bar

# Mature abundance bargraph
ggplot(df %>% 
         filter(var %in% c("sp_N", "C")) %>% 
         mutate(var = factor(var, ordered = TRUE, 
                             levels = c("sp_N", "C"),
                             labels = c("Spawning abundance", "Catch"))),
       aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  scale_fill_manual(values = c("white", "grey70")) +
  theme(legend.position = c(0.15, 0.8),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma, 
                     limits = c(0, max(sp_N$sp_N) * 1.65)) +
  labs(x = NULL, y = "Mature abundance (millions)", linetype = NULL, fill = NULL) -> matn_bar

# Total abundance bargraph
ggplot(df %>% 
         filter(var %in% c("N", "C")) %>% 
         mutate(var = factor(var, ordered = TRUE, 
                             levels = c("N", "C"),
                             labels = c("Mature + immature abundance", "Catch"))),
       aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  scale_fill_manual(values = c("white", "grey70")) +
  theme(legend.position = c(0.2, 0.85),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma, 
                     limits = c(0, max(N$N) * 1.5)) +
  labs(x = NULL, y = "Total abundance (millions)", linetype = NULL, fill = NULL) -> totn_bar

cowplot::plot_grid(matb_bar, matn_bar, totn_bar, align = "hv", nrow = 3) -> bars_her

ggsave(paste0(HERfig_dir, "/biom_abd_barplots.png"), plot = bars_her, dpi = 300, height = 8, width = 7, units = "in")

# Compare LS and HER 

title <- ggdraw() + draw_label("LS")
bars <- cowplot::plot_grid(title, bars, ncol = 1, rel_heights = c(0.1, 1))
title_her <- ggdraw() + draw_label("HER")
bars_her <- cowplot::plot_grid(title_her, bars_her, ncol = 1, rel_heights = c(0.1, 1))

cowplot::plot_grid(bars, bars_her) -> compare_bars
ggsave(paste0(fig_dir, "/compare_biom_abd_barplots.png"), plot = compare_bars, dpi = 300, height = 8, width = 12, units = "in")

# Catch  ----

# For LS:
LS_byyear %>%
  select(year, Catch = tcb, `Spawning biomass` = tot_sp_B_tons, Threshold) %>% 
  full_join(LS_ghl %>% select(year, `Spawning biomass forecast`, GHL, Threshold)) %>% 
  gather("var", "value", -c(year, Threshold)) %>% 
  mutate(var = factor(var, ordered = TRUE, 
                      levels = c("Spawning biomass", "Spawning biomass forecast", "Catch", "GHL"))) -> df

axisb <- tickr(df, year, 5)
ggplot(df, aes(x = year)) +
  geom_bar(aes(y = value, fill = var), colour = "black", size = 0.2, stat = "identity") +
  geom_line(aes(y = Threshold, linetype = "Threshold"), size = 1) +
  scale_fill_manual(values = c("white", "grey85", "grey50", "grey20")) +
  theme(legend.position = c(0.2, 0.8),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma, 
                     breaks = seq(0, 125000, 10000)) +
  labs(x = "", y = "Biomass (tons)\n", linetype = NULL, fill = NULL) -> catch_plot

ggsave(paste0(LSfig_dir, "/biomasscatch_barplot.png"), plot = catch_plot, dpi = 300, height = 4, width = 6, units = "in")

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

ggsave(paste0(fig_dir, "/compare_recruit_plot.png"), plot = recruits, dpi = 300, height = 4, width = 6, units = "in")

# Compare residuals
ggplot(df, aes(x = Year, y = resids)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point() +
  facet_wrap(~ Model, ncol = 1) +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axisr$breaks, labels = axisr$labels) -> rec_resids

ggsave(paste0(fig_dir, "/compare_recruitresids.png"), plot = rec_resids, dpi = 300, height = 5, width = 6, units = "in")

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

ggsave(paste0(fig_dir, "/compare_srcurves.png"), plot = sr_curve, dpi = 300, height = 6, width = 7, units = "in")

# For LS:

# Age-3 recruits compared recruitment estimated from Ricker
LS_byyear %>% 
  ggplot(aes(x = year)) +
  geom_bar(aes(y = init_age_3, colour = "ASA model"),
           stat = "identity", #colour = "grey", 
           fill = "lightgrey",
           width = 0.8, position = position_dodge(width = 0.5)) +
  # If you don't want to show Ricker line, comment out:
  geom_line(aes(y = SR, linetype = "Ricker model"), colour = "black", size = 1) +
  theme(legend.position = c(0.1, 0.75),
        legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.8, 'lines'),
        legend.spacing.y = unit(0, "cm")) +
  # scale_colour_manual(values = "black") +
  scale_colour_manual(values = "grey") +
  scale_x_continuous(breaks = axisb$breaks, labels = axisb$labels) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "\nYear", y = "Age-3 recruits (millions)\n", 
       colour = NULL, linetype = NULL) -> recruits

# residuals
LS_byyear %>% 
  ggplot(aes(x = year, y = res_SR)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = year, xend = year, y = 0, yend = res_SR), 
               size = 0.2, colour = "grey") +
  geom_point() +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) -> resids

# spawner-recruit curve
LS_byyear %>% 
  select(year, SR, init_age_3) %>% 
  filter(year >= min(year) + 3) %>% 
  bind_cols(LS_byyear %>% 
              # shift spawners up by 3 years 
              select(year, tot_sp_B_tons) %>% 
              filter(year <= max(year) - 3) %>% 
              select(- year)) %>% 
  ggplot(aes(x = tot_sp_B_tons)) +
  geom_line(aes(y = SR, linetype = "Ricker model"), 
            colour = "grey", size = 1) +
  geom_point(aes(y = init_age_3, colour = "ASA model")) +
  theme(legend.position = c(0.15, 0.75),
        legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.8, 'lines'),
        legend.spacing.y = unit(0, "cm")) +
  scale_colour_manual(values = "black") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  geom_text_repel(aes(y = init_age_3, label = year), size = 2) +
  labs(x = "\nSpawning stock biomass (tons)", y = "Age-3 recruits (millions)\n",
       linetype = NULL, colour = NULL) -> sr_curve

cowplot::plot_grid(recruits, resids, sr_curve, align = "hv", nrow = 3) -> recruit_plot

ggsave(paste0(LSfig_dir, "/recruit_plot.png"), plot = recruit_plot, dpi = 300, height = 8, width = 6, units = "in")

# For HER:

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

axisr <- tickr(df2, Year, 5)
ggplot(df2, aes(x = Year)) +
  geom_bar(aes(y = age3, colour = "ASA model"),
           stat = "identity", 
           fill = "lightgrey",
           width = 0.8, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = q025, ymax = q975), colour = "black", size = 0.001, width = 0) + 
  # If you don't want to show Ricker line, comment out:
  geom_line(aes(y = SR, linetype = "Ricker model"), colour = "black", size = 1) +
  theme(legend.position = c(0.1, 0.75),
        legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.8, 'lines'),
        legend.spacing.y = unit(0, "cm")) +
  scale_colour_manual(values = "grey") +
  scale_x_continuous(breaks = axisr$breaks, labels = axisr$labels) +
  scale_y_continuous(labels = scales::comma) + 
  theme(legend.position = c(0.15, 0.75)) +
  labs(x = "\nYear", y = "Age-3 recruits (millions)\n", colour = NULL,
       linetype = NULL) -> recruits

# residuals
ggplot(df, aes(x = Year, y = resids)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point() +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axisr$breaks, labels = axisr$labels) -> resids

# spawner-recruit curve
df %>% 
  # Join the spawners to the recruitment years they're associated to (brood year =
  # recruitment year - 3)
  left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                                  D[["mod_nyr"]] + 1, 1),
                       ssb = D[["spawners"]]), by = "Year") %>% 
  # mutate(SR2 = so * ssb * exp(- beta * ssb)) %>% 
  filter(Year >= D[["mod_syr"]] + D[["sage"]]) -> df

# Add in posterior intervals
ricker_sum %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  distinct(Year, q025, q975) %>% 
  left_join(df) -> df

ggplot(df, aes(x = ssb)) +
  geom_ribbon(aes(x = ssb, ymin = q025, ymax = q975),
              alpha = 0.4, fill = "grey80") +
  geom_line(aes(y = SR, linetype = "Ricker model"), 
            colour = "grey", size = 1) +
  geom_point(aes(y = age3, colour = "ASA model")) +
  scale_colour_manual(values = "black") +
  theme(legend.position = c(0.15, 0.75)) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  geom_text_repel(aes(y = age3, label = Year), size = 2) +
  labs(x = "\nSpawning stock biomass (tons)", y = "Age-3 recruits (millions)\n",
       linetype = NULL, colour = NULL) -> sr_curve

cowplot::plot_grid(recruits, resids, sr_curve, align = "hv", nrow = 3, rel_heights=c(1, 0.8, 1)) -> recruit_plot

ggsave(paste0(HERfig_dir, "/recruit_plot.png"), plot = recruit_plot, dpi = 300, height = 8, width = 6, units = "in")

# Egg deposition ----

# Comparison:
df <- as.data.frame(D[["data_egg_dep"]][10:tot_yrs, ])#[seq(tot_yrs - nyr + 1, tot_yrs, 1), ]) 
colnames(df) <- c("Year", "obs", "log_se")
df %>% filter(Year >= D[["mod_syr"]]) %>% 
  # From S. Dressel's bootstrap, have to add new values each year
  mutate(upper = obs + LS_byyear$egg_upper,
         lower = obs - LS_byyear$egg_lower) %>% 
  bind_cols(data.frame(HER = D[["pred_egg_dep"]],
                       LS = LS_byyear$tot_est_egg))%>% 
  gather("Model", "trillions", -c(Year, obs, log_se, upper, lower)) -> df

axisx <- tickr(df, Year, 5)
ggplot(df, aes(x = Year)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "darkgrey", size = 0.001, width = 0.5) + 
  geom_line(aes(y = trillions, colour = Model,
                linetype = Model), size = 1) +  
  geom_point(aes(y = obs, shape = "Historical estimates from survey")) +
  scale_shape_manual(values = 1) +
  scale_color_grey() +
  theme(legend.position = c(0.25, 0.7)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  labs(x = "", y = "Eggs spawned (trillions)\n", shape = "Data") -> eggdep_plot

ggsave(paste0(fig_dir, "/compare_eggdep_plot.png"), plot = eggdep_plot, dpi = 300, height = 4, width = 6, units = "in")

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

ggsave(paste0(fig_dir, "/compare_eggresids.png"), plot = egg_resids, dpi = 300, height = 5, width = 6, units = "in")

# For LS:
LS_byyear %>% 
  rename(Year = year) %>% 
  # 2008 survey egg estimates were extreme high and variable, pmin() helps
  # reduce the scale of the upper confidence interval for better visualization:
  mutate(egg_upper = pmin(tot_obs_egg + egg_upper, max(LS_byyear$tot_obs_egg) * 1.2),
         egg_lower = tot_obs_egg - egg_lower) %>%   ggplot(aes(x = Year)) +
  geom_line(aes(y = tot_est_egg, linetype = "Estimates from model"), size = 1, colour = "grey") +
  scale_colour_manual(values = "grey") +
  geom_point(aes(y = tot_obs_egg, shape = "Historical estimates from survey")) +
  # To get rid of all confidence intervals, comment out following line. To add whiskers, remove width=0
  geom_errorbar(aes(ymin = egg_lower, ymax = egg_upper), colour = "black", width = 0.5, size = 0.001) +
  # theme(legend.position = c(0.25, 0.7)) +
  theme(legend.position = c(0.25, 0.75),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  labs(x = NULL, y = "Eggs spawned (trillions)\n", shape = NULL, linetype = NULL) -> obsfit

# residuals
axis <- tickr(LS_byyear, year, 5)
LS_byyear %>% 
  ggplot(aes(x = year, y = res_tot_egg)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = year, xend = year, y = 0, yend = res_tot_egg), 
               size = 0.2, colour = "grey") +
  geom_point() +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) -> resids

cowplot::plot_grid(obsfit, resids, align = "hv", nrow = 2) -> eggdep_plot

ggsave(paste0(LSfig_dir, "/eggdep_plot.png"), plot = eggdep_plot, dpi = 300, height = 5, width = 6, units = "in")

# For HER with posterior predictive intervals:

# Egg deposition
df <- as.data.frame(D[["data_egg_dep"]])
colnames(df) <- c("Year", "obs_egg", "log_se")
df %>% 
  filter(Year >= D[["mod_syr"]]) %>% 
  mutate(Year = factor(Year)) %>% 
  select(Year, log_se) %>% 
  left_join(egg_sum) -> egg_sum

egg_sum %>%  
  mutate(log_value = log(value),
         # Posterior predictive value
         pp_value = exp(rnorm(n(), mean = log_value, sd = log_se))) -> egg_sum

egg_sum %>% 
  group_by(Year) %>% 
  mutate(pp_mean = mean(pp_value),
         p025 = quantile(pp_value, 0.025),
         p975 = quantile(pp_value, 0.975),
         p250 = quantile(pp_value, 0.250),
         p750 = quantile(pp_value, 0.750)) -> egg_sum

df <- as.data.frame(D[["data_egg_dep"]][ , 1:2])
colnames(df) <- c("Year", "obs")

egg_sum %>% distinct(Year, mean, p025, p975, p250, p750) %>% 
  ungroup() %>% 
  mutate(Year = as.numeric(as.character(Year))) -> egg_sum2

df %>%filter(Year >= D[["mod_syr"]]) -> df

axisx <- tickr(df, Year, 5)

df %>% 
  # 2008 survey egg estimates were extreme high and variable, pmin() helps
  # reduce the scale of the upper confidence interval for better visualization:
  mutate(egg_upper = pmin(obs + LS_byyear$egg_upper, max(obs) * 1.2),
         egg_lower = obs - LS_byyear$egg_lower,
         pred = D[["pred_egg_dep"]]) %>% 
  ggplot(aes(x = Year)) +
  # Posterior predictive intervals 
  geom_point(aes(y = obs, shape = "Historical estimates from survey")) +
  geom_ribbon(data = egg_sum2, aes(x = Year, ymin = p025, ymax = p975),
              alpha = 0.6, fill = "grey70") +
  geom_ribbon(data = egg_sum2, aes(x = Year, ymin = p250, ymax = p750),
              alpha = 0.6, fill = "grey40") +
  geom_line(data = egg_sum2, aes(x = Year, y = mean)) + 
  # Egg dep confidence intervals (To add whiskers, remove width=0)
  geom_errorbar(aes(ymin = egg_lower, ymax = egg_upper), colour = "black", width = 0, size = 0.001) +  scale_colour_manual(values = "grey") +
  theme(legend.position = c(0.25, 0.8),
        legend.spacing.y = unit(0, "cm")) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  scale_y_continuous(limits = c(0, max(df$obs) * 1.2)) +
  labs(x = NULL, y = "Eggs spawned (trillions)\n", shape = NULL, linetype = NULL)  -> obsfit

# residuals
data.frame(Year = D[["year"]],
           resids = D[["resd_egg_dep"]]) %>% 
  ggplot(aes(x = Year, y = resids)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point() +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) -> resids

cowplot::plot_grid(obsfit, resids, align = "hv", nrow = 2, rel_heights = c(1.4, 1)) -> eggdep_plot

ggsave(paste0(HERfig_dir, "/eggdep_plot.png"), plot = eggdep_plot, dpi = 300, height = 5, width = 6, units = "in")

# Survival blocks ----

# Comparison ----

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

survival_blks[c(2:3)] <- lapply(survival_blks[,c(2:3)], prettyNum, digits = 4, big.mark=",")

# Add in LS estimates
df %>%
  select(Year, M = `3`) %>% 
  mutate(survival = exp(-M),
         Model = "HER") %>% 
  bind_rows(LS_byage %>% 
              select(Year, survival) %>% 
              mutate(Model = "LS")) -> df_compare

# Figure summary
tickr(LS_byage, Year, 5) -> axis
ggplot(df_compare, aes(x = Year, y = survival,  colour = Model, linetype = Model)) +
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

ggsave(paste0(fig_dir, "/compare_survival.png"), plot = survival_plot, dpi = 300, height = 4, width = 6, units = "in")

# Plot model comparisons with 95% posterior interval on HER
df_sum %>% 
  ungroup() %>% 
  distinct() %>% 
  melt(id.vars = c("Year", "Blocks"), variable.name = "Model", value.name = "survival") %>%
  left_join(surv_sum %>% distinct(Blocks, q025, q975)) %>% 
  ggplot() +
  annotation_custom(tableGrob(survival_blks, rows = NULL, 
                              theme = ttheme_minimal(base_size = 8, base_colour = "black", base_family = "Times",
                                                     parse = FALSE, #core = list(fg_params = list(hjust = 1)),
                                                     # colhead = list(fg_params = list(hjust = 1)), 
                                                     padding = unit(c(4, 4), "mm"))), 
                    xmin = 1980, xmax = 2015, ymin = 0.1, ymax = 0.35) +
  geom_ribbon(aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.6, fill = "grey70") +
  geom_line(aes(x = Year, y = survival, colour = Model, linetype = Model), size = 1) +
  geom_vline(xintercept = c(1998.5, 2014.5), colour = "lightgrey", linetype = 3, alpha = 0.4) +
  lims(y = c(0, 1)) +
  scale_colour_grey() +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
  labs(x = "", y = "Survival\n") +
  theme(legend.position = c(0.1, 0.8)) -> survival_plot2

ggsave(paste0(fig_dir, "/compare_survival_CI.png"), plot = survival_plot2, dpi = 300, height = 4, width = 6, units = "in")

# For LS:

ggplot(LS_byage, aes(x = Year, y = survival)) +
  geom_vline(xintercept = c(1998.5, 2014.5), colour = "lightgrey", linetype = 3) +
  geom_line(size = 1) +
  geom_point() +
  lims(y = c(0, 1)) +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
  labs(x = "", y = "Survival\n") -> survival_plot

ggsave(paste0(LSfig_dir, "/survival.png"), plot = survival_plot, dpi = 300, height = 4, width = 6, units = "in")

# For HER:

df <- data.frame(D[["year"]], 
                 D[["Mij"]])
colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))

df %>% select(Year, M = `3`) %>% 
  mutate(survival = exp(-M)) -> df

ggplot(df, aes(x = Year, y = survival)) +
  geom_vline(xintercept = c(1998.5, 2014.5), colour = "lightgrey", linetype = 3) +
  geom_line(size = 1) +
  geom_point() +
  lims(y = c(0, 1)) +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
  labs(x = "", y = "Survival\n") -> survival_plot

ggsave(paste0(HERfig_dir, "/survival.png"), plot = survival_plot, dpi = 300, height = 4, width = 6, units = "in")

# HER with 95% and 50% posterior intervals
df_sum %>% 
  ungroup() %>% 
  select(-LS) %>%
  left_join(surv_sum) %>%
  rename(survival = HER) %>%
  distinct(Year, survival, Blocks, q025, q975, q250, q750) -> df_sum2

ggplot(df, aes(x = Year, y = survival)) +
  geom_ribbon(data = df_sum2, aes(x = Year, ymin = q025, ymax = q975),
              alpha = 0.6, fill = "grey70") +
  geom_ribbon(data = df_sum2, aes(x = Year, ymin = q250, ymax = q750),
              alpha = 0.6, fill = "grey40") +
  geom_vline(xintercept = c(1998.5, 2014.5), colour = "lightgrey", linetype = 3) +
  geom_line(size = 1) +
  # geom_point() +
  lims(y = c(0, 1)) +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
  labs(x = "", y = "Survival\n") -> survival_plot_CI

ggsave(paste0(HERfig_dir, "/survival_CI.png"), plot = survival_plot_CI, dpi = 300, height = 4, width = 6, units = "in")

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

ggsave(paste0(fig_dir, "/compare_mat_sel.png"), plot = matsel_plot, dpi = 300, height = 4, width = 6, units = "in")

# For HER with 95% credibility intervals
mat_sum %>% 
  distinct(Blocks, Age, mean, q025, q975) %>% 
  mutate(age = as.numeric(as.character(Age))) %>% 
  ggplot(aes(x = age, y = mean)) + 
  geom_line(aes(linetype = Blocks, group = Blocks)) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  expand_limits(y = 0) +
  labs(x = "\nAge", y = "Proportion\n", linetype = "Time blocks") +
  geom_ribbon(aes(x = age, ymin = q025, ymax = q975, group = Blocks,
                  fill = Blocks),
              alpha = 0.4) +
  scale_fill_manual(values = c("grey70", "grey75"), guide = FALSE) +
  scale_x_continuous(labels = c("3", "4", "5", "6", "7", "8+")) +
  ggtitle("Maturity") +
  theme(legend.position = c(0.7, 0.2),
        plot.title = element_text(hjust = 0.5)) -> mat

# her_matsel %>% filter(param == "Maturity") -> par
# 
# ggplot(par, aes(x = Age, y = proportion)) + 
#   geom_line(aes(linetype = `Time blocks`, group = `Time blocks`)) +
#   geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
#   expand_limits(y = 0) +
#   labs(x = "\nAge", y = "Proportion\n", linetype = "Time blocks") +
#   ggtitle(paste0(par$param[1])) +
#   theme(legend.position = c(0.7, 0.2),
#         plot.title = element_text(hjust = 0.5)) -> mat

sel_sum %>% 
  distinct(Blocks, Age, mean, q025, q975) %>% 
  mutate(Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+"))) %>% 
  # To make sure selectivity is differentiable, it was scaled to have a
  # mean of 1 across all ages. This was done in log space by substracting
  # the mean from the vector of age-specific selectivities. See Tech Doc
  # p 11. Here we normalize it from 0 to 1.
  mutate_if(is.numeric, funs((. - 0) / max(.) - 0)) %>%
  ggplot(aes(x = Age, y = mean)) + 
  geom_line(aes(linetype = Blocks, group = Blocks)) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  expand_limits(y = 0) +
  labs(x = "\nAge", y = NULL, linetype = "Time blocks") +
  geom_ribbon(aes(x = 1:length(Age), ymin = q025, ymax = q975, group = Blocks,
                  fill = Blocks),
              alpha = 0.4) +
  scale_fill_manual(values = c("grey70"), guide = FALSE) +
  ggtitle("Selectivity") +
  theme(legend.position = c(0.7, 0.2),
        plot.title = element_text(hjust = 0.5)) -> sel

cowplot::plot_grid(mat, sel, align = "hv", nrow = 1) -> matsel_plot

ggsave(paste0(HERfig_dir, "/mat_sel.png"), plot = matsel_plot, dpi = 300, height = 4, width = 6, units = "in")

# For LS

ls_matsel %>% filter(param == "Maturity") -> par

ggplot(par, aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, group = `Time blocks`)) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  expand_limits(y = 0) +
  labs(x = "\nAge", y = "Proportion\n", linetype = "Time blocks") +
  ggtitle(paste0(par$param[1])) +
  theme(legend.position = c(0.7, 0.2),
        plot.title = element_text(hjust = 0.5)) -> mat

her_matsel %>% filter(param == "Selectivity") -> par

ggplot(par, aes(x = Age, y = proportion)) + 
  geom_line(aes(linetype = `Time blocks`, group = `Time blocks`)) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 2) +
  expand_limits(y = 0) +
  labs(x = "\nAge", y = NULL, linetype = "Time blocks") +
  scale_y_continuous(breaks = seq(0, max(par$proportion), .25)) +
  ggtitle(paste0(par$param[1])) +
  theme(legend.position = c(0.7, 0.2),
        plot.title = element_text(hjust = 0.5)) -> sel

cowplot::plot_grid(mat, sel, align = "h", nrow = 1) -> matsel_plot

ggsave(paste0(LSfig_dir, "/mat_sel.png"), plot = matsel_plot, dpi = 300, height = 4, width = 6, units = "in")

# Age comps bubbleplot ----

# Bubble plots, just observations (data, no predictions)

# For LS:

axisy <- tickr(LS_byage, Year, 5)
LS_byage %>% 
  melt(id.vars = c("Year", "Age"), measure.vars = c("spawnage_comp_obs", "catchage_comp_obs"), variable.name = "Source", value.name = "proportion") %>% 
  mutate(Source = factor(Source, labels = c("Cast net survey", "Commercial fishery"))) %>% 
  ggplot(aes(x = Age, y = Year, size = proportion)) +
  geom_hline(yintercept = seq(1980, 2010, by = 10), colour = "grey", linetype = 3, alpha = 0.7) +  
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4)) +
  facet_wrap(~ Source) +
  xlab('\nAge') +
  ylab('') +
  guides(size = FALSE) +
  scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) -> agecomps_bubbleplot

ggsave(paste0(LSfig_dir, "/agecomps_bubbleplot.png"), plot = agecomps_bubbleplot, dpi = 300, height = 5, width = 6, units = "in")

# Same figure but using HER:
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

axisy <- tickr(her_agecomps, Year, 5)

ggplot(her_agecomps, aes(x = Age, y = Year, size = obs)) + #*FLAG* could swap size with proportion_scaled
  geom_hline(yintercept = seq(1980, 2010, by = 10), colour = "grey", linetype = 3, alpha = 0.7) +  
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4)) +
  facet_wrap(~ Source) +
  xlab('\nAge') +
  ylab('') +
  guides(size = FALSE) +
  # scale_x_continuous(breaks = unique(df$Age), labels = unique(df$Age)) +
  scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) -> agecomps_bubbleplot

ggsave(paste0(HERfig_dir, "/agecomps_bubbleplot.png"), plot = agecomps_bubbleplot, dpi = 300, height = 5, width = 6, units = "in")

# Residuals age comps ----

# For LS:
LS_byage %>% 
  select(Year, Age, obs = catchage_comp_obs, pred = catchage_comp_est) %>% 
  mutate(Source = "Commercial fishery") %>% 
  bind_rows(LS_byage %>% 
              select(Year, Age, obs = spawnage_comp_obs, pred = spawnage_comp_est) %>% 
              mutate(Source = "Cast net survey")) %>% 
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
         pearson = (obs - pred)/ sqrt(var(pred))) -> df

axis <- tickr(LS_byyear, year, 5)
ggplot(df, aes(x = Age, y = Year, size = resid,
               fill = `Model performance`)) + 
  geom_hline(yintercept = seq(1980, 2010, by = 10), colour = "grey", linetype = 3, alpha = 0.7) +  
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(0, 4)) +
  facet_wrap(~ Source) +
  labs(x = '\nAge', y = '') +
  guides(size = FALSE) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(breaks = axis$breaks, labels = axis$labels) +
  theme(legend.position = "bottom") -> agecomps_residplot

ggsave(paste0(LSfig_dir, "/LS_agecomps_residplot.png"), plot = agecomps_residplot, dpi = 300, height = 5, width = 6, units = "in")

# Same figure for HER

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

ggplot(her_agecomps, aes(x = Age, y = Year, size = resid,
                         fill = `Model performance`)) + 
  geom_hline(yintercept = seq(1980, 2010, by = 10), colour = "grey", linetype = 3, alpha = 0.7) +  
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(0, 3)) +
  facet_wrap(~ Source) +
  labs(x = '\nAge', y = '') +
  guides(size = FALSE) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) +
  theme(legend.position = "bottom") -> agecomps_residplot

ggsave(paste0(HERfig_dir, "/HER_agecomps_residplot.png"), plot = agecomps_residplot, dpi = 300, height = 5, width = 6, units = "in")

# Barplot age comps ----

LS_byage %>% # fishery 
  ggplot() +
  geom_bar(aes(x = Age, y = catchage_comp_obs), 
           stat = "identity", colour = "grey", fill = "lightgrey",
           width = 0.8, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Age, y = catchage_comp_est, group = 1), size = 0.6) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') -> catchage_barplot

ggsave(paste0(LSfig_dir, "/catchage_comps_barplot.png"), plot = catchage_barplot, dpi = 300, height = 8, width = 6, units = "in")

LS_byage %>% # cast net
  ggplot() +
  geom_bar(aes(x = Age, y = spawnage_comp_obs), 
           stat = "identity", colour = "grey", fill = "lightgrey") +
  geom_line(aes(x = Age, y = spawnage_comp_est, group = 1), size = 0.6) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') -> spawnage_barplot

ggsave(paste0(LSfig_dir, "/spawnage_comps_barplot.png"), plot = spawnage_barplot, dpi = 300, height = 8, width = 6, units = "in")

# Same figs for HER
her_agecomps %>% # fishery
  filter(Source == "Commercial fishery") %>% 
  ggplot() + 
  geom_bar(aes(x = Age, y = obs), 
           stat = "identity", colour = "grey", fill = "lightgrey",
           width = 0.8, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Age, y = pred, group = 1), size = 0.6) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Commercial fishery age compositions") -> catchage_barplot

ggsave(paste0(HERfig_dir, "/catchage_comps_barplot.png"), plot = catchage_barplot, dpi = 300, height = 8, width = 6, units = "in")

# Barplot with 95% posterior predictive intervals
cm_comp_ppi <- cm_comp_ppi %>% 
  mutate(age = as.numeric(as.character(Age)),
         Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+")),
         Cohort = factor(ifelse(age < 8, Year - age, "plus"), ordered = TRUE))

her_agecomps %>% # fishery
  filter(Source == "Commercial fishery") %>% 
  left_join(cm_comp_ppi) %>% #View()
  ggplot() + 
  geom_bar(aes(x = Age, y = obs), 
           stat = "identity", colour = "grey", fill = "lightgrey",
           width = 0.8, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = Age, y = pred), size = 1) +
  geom_errorbar(aes(x = Age, ymin = q025, ymax = q975),
                colour = "black", width = 0, size = 0.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Commercial fishery age compositions") -> catchage_barplot2

ggsave(paste0(HERfig_dir, "/catchage_comps_barplotPPI.png"), plot = catchage_barplot2, dpi = 300, height = 8, width = 6, units = "in")

# colour palette for tracking cohorts
mycols <- c("#f4cc70","#de7a22", "#20948b", "#6ab187", "#9a9eab", "white")#BCBABE")
nages <- D[["nage"]] - D[["sage"]]
xtra <- nages * ((n_distinct(cm_comp_ppi$Cohort)-1)/nages - floor((n_distinct(cm_comp_ppi$Cohort)-1)/nages))
cohort_fill <- c(rep(mycols[1:5], (n_distinct(cm_comp_ppi$Cohort)-1)/nages), mycols[1:xtra], mycols[6])
cohort_cols <- replace(cohort_fill, cohort_fill == "white", "black") # Plus group always white bar with black border

her_agecomps %>% # fishery
  filter(Source == "Commercial fishery") %>% 
  left_join(cm_comp_ppi) %>%
  ggplot() + 
  geom_bar(aes(x = Age, y = obs, colour = Cohort, fill = Cohort), 
           stat = "identity", show.legend = FALSE,
           width = 0.8, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = cohort_fill) +
  scale_colour_manual(values = cohort_cols) +
  geom_point(aes(x = Age, y = pred), size = 1) +
  geom_errorbar(aes(x = Age, ymin = q025, ymax = q975),
                colour = "black", width = 0, size = 0.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Commercial fishery age compositions") -> catchage_barplot3

ggsave(paste0(HERfig_dir, "/catchage_comps_colourbarplotPPI.png"), plot = catchage_barplot3, dpi = 300, height = 8, width = 7, units = "in")

# survey
her_agecomps %>% 
  filter(Source == "Cast net survey") %>% 
  ggplot() + 
  geom_bar(aes(x = Age, y = obs), 
           stat = "identity", colour = "grey", fill = "lightgrey",
           width = 0.8, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Age, y = pred, group = 1), size = 0.6) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Cast net survey age compositions") -> spawnage_barplot

ggsave(paste0(HERfig_dir, "/spawnage_comps_barplot.png"), plot = spawnage_barplot, dpi = 300, height = 8, width = 6, units = "in")

# Barplot with 95% posterior predictive intervals
sp_comp_ppi <- sp_comp_ppi %>% 
  mutate(age = as.numeric(as.character(Age)),
         Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+")),
         Cohort = factor(ifelse(age < 8, Year - age, "plus"), ordered = TRUE))

her_agecomps %>% # survey
  filter(Source == "Cast net survey") %>% 
  left_join(sp_comp_ppi) %>% #View()
  ggplot() + 
  geom_bar(aes(x = Age, y = obs), 
           stat = "identity", colour = "grey", fill = "lightgrey",
           width = 0.8, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = Age, y = pred), size = 1) +
  geom_errorbar(aes(x = Age, ymin = q025, ymax = q975),
                colour = "black", width = 0, size = 0.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') + 
  ggtitle("Cast net survey age compositions") -> spawnage_barplot2

ggsave(paste0(HERfig_dir, "/spawnage_comps_barplotPPI.png"), plot = spawnage_barplot2, dpi = 300, height = 8, width = 6, units = "in")

her_agecomps %>% # survey
  filter(Source == "Cast net survey") %>% 
  left_join(sp_comp_ppi) %>%
  ggplot() + 
  geom_bar(aes(x = Age, y = obs, colour = Cohort, fill = Cohort), 
           stat = "identity", show.legend = FALSE,
           width = 0.8, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = cohort_fill) +
  scale_colour_manual(values = cohort_cols) +
  geom_point(aes(x = Age, y = pred), size = 1) +
  geom_errorbar(aes(x = Age, ymin = q025, ymax = q975),
                colour = "black", width = 0, size = 0.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') +
  ggtitle("Cast net survey age compositions") -> spawnage_barplot3

ggsave(paste0(HERfig_dir, "/spawnage_comps_colourbarplotPPI.png"), plot = spawnage_barplot3, dpi = 300, height = 8, width = 6, units = "in")

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
  theme(legend.position = "bottom") + 
  ggtitle("Cast net survey age compositions")  -> spcomp_barplot

ggsave(paste0(fig_dir, "/compare_spcomp_barplot.png"), plot = spcomp_barplot, dpi = 300, height = 8, width = 6, units = "in")

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
  theme(legend.position = "bottom") + 
  ggtitle("Commercial fishery age compositions") -> catchcomp_barplot

ggsave(paste0(fig_dir, "/compare_catchcomp_barplot.png"), plot = catchcomp_barplot, dpi = 300, height = 8, width = 6, units = "in")

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

df <- as.data.frame(D[["data_sp_waa"]])
colnames(df) <- c("Year", paste(D[['sage']] : D[['nage']]))
df %>% 
  gather("Age", "weight", -Year) %>% 
  filter(Year >= D[["mod_syr"]]) %>% 
  group_by(Age) %>% 
  summarise(waa = mean(weight)) -> waa_bar

df <- as.data.frame(D[["mat"]])
colnames(df) <- c(paste(D[['sage']] : D[['nage']]))
df %>% 
  mutate(Year = D[["year"]]) %>% 
  gather("Age", "maturity", -Year) %>% 
  filter(Year >= D[["mod_syr"]]) %>% 
  group_by(Age) %>% 
  summarise(mat = mean(maturity)) -> mat_bar

df <- as.data.frame(D[["Mij"]])
colnames(df) <- c(paste(D[['sage']] : D[['nage']]))
df %>% 
  mutate(Year = D[["year"]]) %>% 
  gather("Age", "mort", -Year) %>% 
  filter(Year >= D[["mod_syr"]]) %>% 
  #group_by(Age) %>% 
  summarise(mort = mean(mort)) %>% 
  pull(mort) -> mort_bar

lx <- c(1:6)

lx[1] <- 1

for(i in 2:6){
  lx[i] <- lx[i-1] * exp(-mort_bar)
}
# plus group
lx[6] = lx[6] / (1 - exp(-mort_bar))
lx  

waa_bar %>% 
  left_join(mat_bar) %>% 
  mutate(lx = lx,
         fec = waa * mat,
         phi = lx * fec) %>% 
  summarize(phie = sum(phi)) %>% 
  pull(phie) -> phie

# From .par file 
log_ro <- D[["theta"]][4] # theta(4)
log_reck <- D[["theta"]][5]  # theta(5)

ro <- exp(log_ro) # equilibrium recruitment
reck <-  exp(log_reck) + 1.0 # recruitment compensation ratio

# stock-recruit parameters
so <- reck/phie;
beta <- log(reck) / (ro * phie) # beta = (reck - 1.0)/(ro *phiE)

so; beta
D[["so"]];D[["beta"]]

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
