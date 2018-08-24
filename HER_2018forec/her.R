
# HER - ADFG's new herring model. Original ADMB code written by SJD Martell. Helper
# files and documentation contributed by M. Rudd.

# Libraries/source files ----

library(R2admb) # run ADMB in R
library(tidyverse)
library(data.table)

source("R/tools.r")
source("R/helper.r") 

# LS results ----

# LS 2018 forecast results
LS_forec <- read_csv("HER_2018forec/LS_2018forec_results.csv")
LS_byage <- read_csv("HER_2018forec/LS_2018forec_results_byage.csv")
LS_byyear <- read_csv("HER_2018forec/LS_2018forec_results_byyear.csv")

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
LS_yearminus <- read_csv("HER_2018forec/yearminus_matbio.csv")

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

  # Run model in ADMB ----

# Can't figure out how to put file path directly from project root into the
# compile_admb() function
#setwd(file.path(getwd(), "HER_2018forec"))
setwd(file.path(getwd(), "HER"))

# Running model in R
setup_admb()
compile_admb("her", verbose = TRUE)
run_admb("her", verbose = TRUE) 

# Running her model in command line
# 1) Open command prompt ("C:/ADMB/ADMB Command Prompt (MinGW 64Bit)")
# 2) Navigate to AlaskaHerring deep inside S drive using 
# > cd ../..
# > S:
# > cd "S:\Region1Shared-DCF\Research\Herring-Dive Fisheries\Herring\ADMB Rudd Workshop 2018\AlaskaHerring\HER"
# Compile
# > admb her 
# Run
# > her

# Diagnostics ----
P <- read_fit("her")
P[["nopar"]]
P[["nlogl"]]
P[["logDetHess"]]
P[["maxgrad"]]

# Results ----

D <- read_admb("her")

D$fore_sb
D$ghl


# Spawning biomass ----

# post-fishery

tot_yrs <- D[["dat_nyr"]] - D[["dat_syr"]] + 1


df <- data.frame(Year = D[["year"]], 
                 spB = D[["ssb"]] / 0.90718, # convert to short tons
                 catch = D[["data_catch"]][10:47, 2]#[nyr + 1, tot_yrs, 1), 2] # just the column of catch, already in short tons
) %>% 
  mutate(matB = spB + catch,
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
  # geom_point(data = df, aes(x = Year, y = tons, colour = Model, linetype = Model, shape = Model), size = 1) +
  geom_line(data = df, aes(x = Year, y = tons, colour = Model, linetype = Model), size = 1) +
  geom_point(data = srv_index, aes(x = Year, y = srv_spB, shape = "Historical estimates from survey")) +
  scale_shape_manual(values = 1) +
  scale_colour_grey() +
  theme(legend.position = c(0.25, 0.7)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "", y = "Spawning biomass (tons)\n", shape = "Data") -> spbiomass_plot

ggsave("figs/compare_spbiomass_plot.png", plot = spbiomass_plot, dpi = 300, height = 4, width = 6, units = "in")

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
tickryr <- data.frame(Year = 1980:2020)
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
  scale_y_continuous(limits = c(0, max(LS_yearminus$matbio_tons) * 1.1),
                     labels = scales::comma) +
  labs(x = "", y = "Spawning biomass (tons)\n", shape = NULL) -> past_matbio

ggsave("figs/LS/compare_past_matbio.png", plot = past_matbio, dpi = 300, height = 4, width = 6, units = "in")

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

ggsave("figs/LS/biomasscatch_barplot.png", plot = catch_plot, dpi = 300, height = 4, width = 6, units = "in")

# For HER:

# Wait on this b/c currently mature biomass is not forecasted.

# Forecasted comps ----

# Comparison:

# For LS:

LS_forec %>% 
  mutate(biom = prettyNum(for_mat_baa_tons, digits = 1, big.mark=","),
         `Weight-at-age (*g*)` = formatC(round(for_waa, 1), small.interval = 1),
         `Proportion mature` = formatC(round(for_mat_prop, 2), small.interval = 2),
         biom = ifelse(for_mat_baa_tons == max(for_mat_baa_tons), 
                                         paste0("**", biom, "**"), biom)) %>% 
  select(Age, `Mature biomass (*t*)` = biom, `Weight-at-age (*g*)`, `Proportion mature`) -> forec_age

write_csv
# For HER:

# Wait on this b/c HER doesn't currently output forecast numbers or biomass at
# age

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

ggsave("figs/compare_recruit_plot.png", plot = recruits, dpi = 300, height = 4, width = 6, units = "in")

# Compare residuals
ggplot(df, aes(x = Year, y = resids)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point() +
  facet_wrap(~ Model, ncol = 1) +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axisr$breaks, labels = axisr$labels) -> rec_resids

ggsave("figs/compare_recruitresids.png", plot = rec_resids, dpi = 300, height = 5, width = 6, units = "in")

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

ggsave("figs/compare_srcurves.png", plot = sr_curve, dpi = 300, height = 6, width = 7, units = "in")

# For LS:

# Age-3 recruits compared recruitment estimated from Ricker
LS_byyear %>% 
  ggplot(aes(x = year)) +
  geom_line(aes(y = SR, linetype = "Ricker model"), colour = "grey", size = 1) +
  geom_point(aes(y = init_age_3, colour = "ASA model")) +
  theme(legend.position = c(0.15, 0.75),
        legend.key = element_rect(size = 0.5),
        legend.key.size = unit(0.8, 'lines'),
        legend.spacing.y = unit(0, "cm")) +
  scale_colour_manual(values = "black") +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
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

ggsave("figs/LS/recruit_plot.png", plot = recruit_plot, dpi = 300, height = 8, width = 6, units = "in")

# For HER:

# Age-3 recruits compared recruitment estimated from Ricker
df <- data.frame(Year = D[["years"]],
                 age3 = D[["Nij"]][, 1]) %>% 
  left_join(data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                            D[["mod_nyr"]] + 1, 1),
                 SR = D[["recruits"]],
                 resids = D[["resd_rec"]]), by = "Year")

axisr <- tickr(df, Year, 5)
ggplot(df, aes(x = Year)) +
  geom_line(aes(y = SR, linetype = "Ricker model"), colour = "grey", size = 1) +
  geom_point(aes(y = age3, colour = "ASA model")) +
  scale_colour_manual(values = "black") +
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
  mutate(SR2 = so * ssb * exp(- beta * ssb)) %>% 
  filter(Year >= D[["mod_syr"]] + D[["sage"]]) -> df

ggplot(df, aes(x = ssb)) +
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

cowplot::plot_grid(recruits, resids, sr_curve, align = "hv", nrow = 3) -> recruit_plot

ggsave("figs/HER/recruit_plot.png", plot = recruit_plot, dpi = 300, height = 8, width = 6, units = "in")

# Egg deposition ----

# Comparison:
df <- as.data.frame(D[["data_egg_dep"]][10:47, ])#[seq(tot_yrs - nyr + 1, tot_yrs, 1), ]) 
colnames(df) <- c("Year", "obs", "log_se")
df %>% filter(Year >= D[["mod_syr"]]) %>% 
  bind_cols(data.frame(HER = D[["pred_egg_dep"]],
                       LS = LS_byyear$tot_est_egg[1:38]))%>% 
  gather("Model", "trillions", -c(Year, obs, log_se)) -> df

ggplot(df, aes(x = Year)) +
  geom_line(aes(y = trillions, colour = Model,
                linetype = Model), size = 1) +  
  geom_point(aes(y = obs, shape = "Historical estimates from survey")) +
  scale_shape_manual(values = 1) +
  scale_color_grey() +
  theme(legend.position = c(0.25, 0.7)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  labs(x = "", y = "Eggs spawned (trillions)\n", shape = "Data") -> eggdep_plot

ggsave("figs/compare_eggdep_plot.png", plot = eggdep_plot, dpi = 300, height = 4, width = 6, units = "in")

# Compare resids

data.frame(Year = D[["year"]],
                 resids = D[["resd_egg_dep"]],
                 Model = "HER") %>% 
  bind_rows(LS_byyear %>% 
              select(Year = year, resids = res_tot_egg) %>% 
              mutate(Model = "LS")) %>% 
  ggplot(aes(x = Year, y = resids)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point() +
  facet_wrap(~ Model, ncol = 1) +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) -> egg_resids

ggsave("figs/compare_eggresids.png", plot = egg_resids, dpi = 300, height = 5, width = 6, units = "in")

# For LS:
LS_byyear %>% 
  rename(Year = year) %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = tot_est_egg, linetype = "Estimated from model"), size = 1, colour = "grey") +
  scale_colour_manual(values = "grey") +
  geom_point(aes(y = tot_obs_egg, shape = "Historical estimates from survey")) +
  theme(legend.position = c(0.25, 0.7)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  labs(x = NULL, y = "Eggs spawned (trillions)\n", shape = NULL, linetype = NULL) -> obsfit

# residuals
LS_byyear %>% 
  ggplot(aes(x = year, y = res_tot_egg)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = year, xend = year, y = 0, yend = res_tot_egg), 
               size = 0.2, colour = "grey") +
  geom_point() +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) -> resids

cowplot::plot_grid(obsfit, resids, align = "hv", nrow = 2) -> eggdep_plot

ggsave("figs/LS/eggdep_plot.png", plot = eggdep_plot, dpi = 300, height = 5, width = 6, units = "in")

# For HER:

df <- as.data.frame(D[["data_egg_dep"]][ , 1:2])
colnames(df) <- c("Year", "obs")
df %>% 
  filter(Year >= D[["mod_syr"]]) %>% 
  mutate(pred = D[["pred_egg_dep"]]) %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = pred, linetype = "Estimated from model"), size = 1, colour = "grey") +
  scale_colour_manual(values = "grey") +
  geom_point(aes(y = obs, shape = "Historical estimates from survey")) +
  theme(legend.position = c(0.25, 0.7)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  labs(x = NULL, y = "Eggs spawned (trillions)\n", shape = NULL, linetype = NULL) -> obsfit

# residuals
data.frame(Year = D[["year"]],
                 resids = D[["resd_egg_dep"]]) %>% 
  ggplot(aes(x = Year, y = resids)) + 
  geom_hline(yintercept = 0, colour = "grey", size = 1) +
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = resids), 
               size = 0.2, colour = "grey") +
  geom_point() +
  labs(x = "\nYear", y = "Residuals\n") +
  scale_x_continuous(breaks = axis$breaks, labels = axis$labels) -> resids

cowplot::plot_grid(obsfit, resids, align = "hv", nrow = 2) -> eggdep_plot

ggsave("figs/HER/eggdep_plot.png", plot = eggdep_plot, dpi = 300, height = 5, width = 6, units = "in")

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

ggsave("figs/compare_survival.png", plot = survival_plot, dpi = 300, height = 4, width = 6, units = "in")

# Maturity ----

df <- data.frame(D[["year"]], 
                 D[["mat"]])
colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
write.table(df,sep="\t", row.names=FALSE)

df %>% 
  gather("Age", "maturity", -Year) %>% 
  mutate(Model = "HER",
         Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+"))) %>% 
  bind_rows(LS_byage %>% 
              select(Year, Age, maturity) %>% 
              mutate(Model = "LS")) %>% 
  group_by(Model, Age, maturity) %>% 
  summarise(min = min(Year),
            max = max(Year)) %>% 
  mutate(`Maturity blocks` = paste0(min, "-", max),
         combos = paste0(Model, " ", `Maturity blocks`)) -> df

ggplot(df, aes(x = Age, y = maturity, colour = Model)) + 
  geom_line(aes(linetype = `Maturity blocks`,  group = combos)) +
  scale_colour_grey() +
  lims(y = c(0, 1)) +
  labs(x = "\nAge", y = "Maturity\n") +
  theme(legend.position = c(0.8, 0.6)) -> maturity_plot

ggsave("figs/compare_maturity.png", plot = maturity_plot, dpi = 300, height = 4, width = 6, units = "in")

# Selectivity ----

# number of yrs in the model 
nyr <- D[["mod_nyr"]] - D[["mod_syr"]] + 1

df <- data.frame(D[["year"]], 
                 D[["Sij"]][1:nyr, ]) # change number of rows
colnames(df) <- c("Year", paste(D[['sage']]:D[['nage']]))
# write.table(df,sep="\t", row.names=FALSE)

df %>% 
  gather("Age", "selectivity", -Year) %>% 
  mutate(Model = "HER",
         Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+")),
         # To make sure selectivity is differentiable, it was scaled to have a
         # mean of 1 across all ages. This was done in log space by substracting
         # the mean from the vector of age-specific selectivities. See Tech Doc
         # p 11. *FLAG* here's an attempt to normalize it from 0 to 1.
         selectivity2 = ( selectivity - min(selectivity) ) / ( max(selectivity) - min(selectivity) ) ) %>%
  bind_rows(LS_byage %>% 
              select(Year, Age, selectivity = gear_select) %>% 
              mutate(Model = "LS")) %>% 
  group_by(Model, Age, selectivity, selectivity2, selectivity3) %>% 
  summarise(min = min(Year),
            max = max(Year)) %>% 
  mutate(`Selectivity blocks` = paste0(min, "-", max),
         combos = paste0(Model, " ", `Selectivity blocks`)) -> df

ggplot(df, aes(x = Age, y = selectivity, colour = Model)) + 
  geom_line(aes(linetype = `Selectivity blocks`, group = combos)) +
  geom_line(aes(x = Age, y = selectivity2, group = Model), colour = "black", linetype = 2) +
  scale_colour_grey() +
  lims(y = c(0, 1.8)) +
  labs(x = "\nAge", y = "Selectivity\n") +
  theme(legend.position = c(0.8, 0.25)) -> selectivity_plot

ggsave("figs/compare_selectivity.png", plot = selectivity_plot, dpi = 300, height = 4, width = 6, units = "in")

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

ggsave("figs/LS/catchcomps_bubbleplot.png", plot = agecomps_bubbleplot, dpi = 300, height = 5, width = 6, units = "in")

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

ggsave("figs/HER/catchcomps_bubbleplot.png", plot = agecomps_bubbleplot, dpi = 300, height = 5, width = 6, units = "in")

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
         `Model performance` = ifelse(raw >= 0, "Predicted less than observed", "Predicted greater than observed"),
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

ggplot(df, aes(x = Age, y = Year, size = pearson,
             fill = `Model performance`)) + 
  geom_hline(yintercept = seq(1980, 2010, by = 10), colour = "grey", linetype = 3, alpha = 0.7) +  
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(0, 4)) +
  facet_wrap(~ Source) +
  labs(x = '\nAge', y = '') +
  guides(size = FALSE) +
  scale_fill_manual(values = c("white", "black")) +
  scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) +
  theme(legend.position = "bottom") -> agecomps_residplot

ggsave("figs/LS/agecomps_residplot.png", plot = agecomps_residplot, dpi = 300, height = 5, width = 6, units = "in")

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
         `Model performance` = ifelse(raw >= 0, "Predicted less than observed", "Predicted greater than observed"),
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

ggplot(her_agecomps, aes(x = Age, y = Year, size = pearson,
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

ggsave("figs/HER/agecomps_residplot.png", plot = agecomps_residplot, dpi = 300, height = 5, width = 6, units = "in")

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

ggsave("figs/LS/catchage_comps_barplot.png", plot = catchage_barplot, dpi = 300, height = 8, width = 6, units = "in")

LS_byage %>% # cast net
  ggplot() +
  geom_bar(aes(x = Age, y = spawnage_comp_obs), 
           stat = "identity", colour = "grey", fill = "lightgrey") +
  geom_line(aes(x = Age, y = spawnage_comp_est, group = 1), size = 0.6) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') -> spawnage_barplot

ggsave("figs/LS/spawnage_comps_barplot.png", plot = spawnage_barplot, dpi = 300, height = 8, width = 6, units = "in")

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
  labs(x = '\nAge', y = 'Proportion-at-age\n') -> catchage_barplot

ggsave("figs/HER/catchage_comps_barplot.png", plot = catchage_barplot, dpi = 300, height = 8, width = 6, units = "in")

her_agecomps %>% # survey
  filter(Source == "Cast net survey") %>% 
  ggplot() + 
  geom_bar(aes(x = Age, y = obs), 
           stat = "identity", colour = "grey", fill = "lightgrey",
           width = 0.8, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Age, y = pred, group = 1), size = 0.6) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') -> spawnage_barplot

ggsave("figs/HER/spawnage_comps_barplot.png", plot = spawnage_barplot, dpi = 300, height = 8, width = 6, units = "in")

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

ggsave("figs/compare_spcomp_barplot.png", plot = spcomp_barplot, dpi = 300, height = 8, width = 6, units = "in")


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

ggsave("figs/compare_catchcomp_barplot.png", plot = catchcomp_barplot, dpi = 300, height = 8, width = 6, units = "in")

# Steve's stuff ----

# Read in the data from the model report, par, and cor files.
# D <- read.admb("HER_2018forec/her") # read.admb() from globals.r
D <- read_admb("her") # read.admb() from globals.r

# sb.file <- "HER_2018forec/ssb.ps" # Will only exist if you've run -mceval
sb.file <- "ssb.ps" # Will only exist if you've run -mceval
if(file.exists(sb.file)){
  D$post.samp.ssb=read.table(sb.file)
  colnames(D$post.samp.ssb) <- paste0("year",D$year)
}


plot.catch <- function(D = D, nm = "data_ct_raw",...) {
  df <- as.data.frame(D[[nm]])
  colnames(df) <- c("Year","Catch","log.se")
  z  <- 1.96
   df <- df %>% 
    mutate(ln.ct = log(Catch),
           lci = exp(ln.ct - z*log.se),
           uci = exp(ln.ct + z*log.se),
           std = 1.96 * sqrt(log(log.se+1)),
           lower = Catch - std * Catch,
           upper = Catch + std * Catch)
  
  tickr(df, Year, 5)
  ggplot(df,aes(Year,Catch)) +
    geom_point() +
    geom_pointrange(aes(ymin = lci, ymax = uci),size=0.5,fatten=2) + 
    labs(x="Year")
}


plot.waa <- function(D=D, nm = "data_sp_waa",...) {
  df <- as.data.frame(D[[nm]])
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  gdf <- 	gather(df,"Age","Weight.at.age",-Year) %>% 
    transform(Age=as.integer(Age)) %>% 
    mutate(Cohort=as.factor(Year-Age))
  
  ggplot(gdf,aes(Year,Weight.at.age,color=Cohort)) + 
    geom_line(alpha=0.90) + 
    geom_point(alpha=0.5,aes(fill=Cohort),show.legend=FALSE,size=0.5) +
    labs(x="Year",y="Weight-at-age (grams)",color="Cohort") +
    guides(col = guide_legend(ncol = 9)) +
    theme(legend.position="bottom") +ggtitle(D$Model)
}

plot.comp <- function(D=D, nm = "data_cm_comp",...) {
  df <- as.data.frame(D[[nm]])
  df[df==-9] <- NA
  colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
  gdf <- 	gather(df,"Age","Proportion",-Year) %>%
    transform(Age=as.integer(Age)) %>%
    mutate(Cohort=as.factor(Year-Age))
  
  ggplot(gdf,aes(Year,Age,color=Cohort)) + 
    geom_point(alpha=0.5,aes(size=abs(Proportion)),show.legend=FALSE) +
    scale_size_area(max_size=8) + 
    labs(x="Year",...) + ggtitle(D$Model)
}

plot.ssb <- function(D=D){
  #qplot(D$year,D$ssb/1000,geom="line") + ylim(c(0,NA)) +
  #labs(x="Year",y="Female Spawning Stock Biomass (1000 mt)")
  
  # df <- data.frame(year=seq(D$mod_syr,D$mod_nyr),ssb=D$ssb)
  # ggplot(df,aes(year,ssb/1000)) +geom_line() + ylim(c(0,NA)) + 
  # labs(x="Year",y="Female Spawning Stock Biomass (1000 mt)") +
  # ggtitle(D$Model)
  
  
  id <- grep("sd_ssb",D$fit$names)
  
  ssb.df <- data.frame(year=seq(D$mod_syr,D$mod_nyr),
                       SSB = D$fit$est[id]/1000,
                       sdSSB = D$fit$std[id]/1000) %>% 
    mutate(lci=SSB-1.96*sdSSB,uci=SSB+1.96*sdSSB)
  
  ggplot(ssb.df,aes(year,SSB)) + 
    geom_line() +
    geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.15)+
    labs(x="Year",y="Female Spawning Stock Biomass (1000 mt)") +
    ggtitle(D$Model)
  
  
}

plot.resd <- function(D=D, nm = "resd_cm_comp", ...) {
  
  # Dealing with composition data.
  if( grepl("comp",nm) ){
    df <- as.data.frame(cbind(D[['mod_syr']]:D[['mod_nyr']],D[[nm]]))
    # df[df==-9] <- NA
    colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
    gdf <- 	gather(df,"Age","Residual",-Year) %>%
      transform(Age=as.integer(Age)) %>%
      mutate(Cohort=as.factor(Year-Age))
    
    ggplot(gdf,aes(Year,Age,color=factor(sign(Residual)))) + 
      geom_point(alpha=0.5,aes(size=abs(Residual)),show.legend=TRUE) +
      scale_size_area(max_size=8) + 
      labs(x="Year",color="Sign",size="Residual",...)+
      ggtitle(D$Model)
    
  } else if( grepl("egg",nm) ) {
    df <- as.data.frame(cbind(D[['year']],D[[nm]]))
    colnames(df) <- c("Year","Residual")
    
    ggplot(df,aes(Year,Residual)) + geom_point() +
      geom_segment(aes(x = Year, xend = Year, y = 0, 
                       yend = Residual),data=df,size=0.2) +
      labs(y="Residual (egg deposition)")+
      ggtitle(D$Model)
    
  } else if( grepl("rec",nm) ) {
    df <- as.data.frame(cbind(D[['rec_years']],D[[nm]]))
    colnames(df) <- c("Year","Residual")
    
    ggplot(df,aes(Year,Residual)) + geom_point() +
      geom_segment(aes(x = Year, xend = Year, y = 0, 
                       yend = Residual),data=df,size=0.2)+
      labs(y="Residual (recruitment deviation)")+
      ggtitle(D$Model)
  } else if( grepl("catch",nm) ) {
    df <- as.data.frame(cbind(D[['year']],D[[nm]]))
    colnames(df) <- c("Year","Residual")
    
    ggplot(df,aes(Year,Residual)) + geom_point() +
      geom_segment(aes(x = Year, xend = Year, y = 0, 
                       yend = Residual),data=df,size=0.2)+
      labs(y="Residual (commercial catch)") +
      ggtitle(D$Model)
  }
}

plot.ft <- function(D) {
  id <- grep("log_ft_pars",D$fit$names)
  log.ft.mle <- D$fit$est[id]
  log.ft.std <- D$fit$std[id]
  
  df <- data.frame(Year=D$year,Ft = exp(log.ft.mle),
                   lci = exp(log.ft.mle-1.96*log.ft.std),
                   uci = exp(log.ft.mle+1.96*log.ft.std))
  
  ggplot(df,aes(Year,Ft)) + geom_line() + 
    geom_ribbon(aes(x=Year,ymin=lci,ymax=uci),alpha=0.15)+
    labs(y="Instantaneous fishing mortality") + 
    ggtitle(D$Model)
  
}

plot.ft.post <- function(D=D) {
  
  # if(is.null(D$post.samp)) return()
  ps <- D$post.samp
  ix <- sample(1:ncol(ps),1000,replace=TRUE)
  ps <- ps[ix,]
  colnames(ps) <- D$fit$names[1:ncol(ps)]	
  
  # select the log_ft_pars columns
  px <- ps[,grepl("log_ft_pars",colnames(ps))]
  yr <- seq(D$mod_syr,D$mod_nyr)	
  px <- as.data.frame(px)
  colnames(px) <- paste(yr)
  
  # gather
  gx <- gather(px,Year,Value)
  # plot
  ggplot(gx,aes(Year,exp(Value))) +
    geom_violin(alpha=0.25,fill="red",size=0.15,
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    labs(x="Year",y="Average fishing mortality rate (ft)")+
    ylim(c(0,NA))+
    scale_x_discrete(breaks=seq(D$mod_syr,D$mod_nyr,by=5))
  
}

# Spawning stock biomass posterior samples
plot.ssb.post <- function(D = D) {
  ssb.ps <- D$post.samp.ssb
  colnames(ssb.ps) <- paste(D$year)
  df <- gather(ssb.ps, Year, SSB) %>% 
    group_by(Year) %>% 
    summarize(median = median(SSB) / 1000,
              q025 = quantile(SSB, 0.025) / 1000,
              q975 = quantile(SSB, 0.975) / 1000)
  
  # ggplot(df, aes(Year, SSB / 1000)) + 
  #   geom_violin(alpha = 0.25, fill = "steel blue", size = 0.25,
  #               draw_quantiles = c(0.25, 0.5, 0.75)) +
  #   ylim(c(0, NA)) + 
  #   labs(x = "Year", y = "Spawning Stock Biomass (1000 t)\n") +
  #   scale_x_discrete(breaks = seq(D$mod_syr, D$mod_nyr, by = 5))

  ggplot(df, aes(x = Year, y = SSB / 1000)) + 
    geom_line() +
    ylim(c(0, NA)) + 
    labs(x = "Year", y = "Spawning Stock Biomass (1000 t)\n") +
    scale_x_discrete(breaks = seq(D$mod_syr, D$mod_nyr, by = 5))
  
  ggplot() +
    # geom_point( aes(x = Year, y = median)) +
    # geom_line( aes(x = Year, y = median, group = 1)) +
    geom_ribbon(data = df, 
                aes(x = Year, ymin = q025, ymax = q975),
                colour = "black", fill = "grey") #+, 
                # alpha = 0.2, fill = "grey") +
    # labs(x = "\nYear",
    #      y = "Spawning stock biomass (1000 t)\n") +
    # # theme(legend.position = "none") +
    # xlim(c(0, 130)) 
    
  # }
}

# Egg deposition 
plot.eggdepfit <- function(D = D, sfx = "egg_dep", fit = TRUE) {
  data <- paste0("data_", sfx)
  data <- as.data.frame(D[[data]])
  colnames(data) <- c("year","index","log.se")
  
  pred <- paste0("pred_",sfx)
  pred <- as.data.frame(cbind(D$year,D[[pred]]))
  colnames(pred) <- c("year","pred")
  
  # From Sherri's bootstrap, have to add new values each year
  U <- c(1.18, 1.12, 1.10, 0.95, 1.23, 1.13, 0.93, 1.67, 2.51, 0.98, 0.80, 1.259686251, 1.851147636, 1.751126689, 0.560987576, 1.508697833, 1.633749193, 1.695401525, 1.255367509, 2.387620843, 2.795968238, 2.215761696, 1.462234716, 2.289501604, 2.650921062, 6.384923885, 2.279693146, 2.872760889, 29.05686308, 3.863145759, 4.816134565, 4.222205391, 1.634805164, 4.043867944, 1.288746439, 1.332721825, 2.445122264,1.69865191)
  L <- c(1.18, 1.12, 1.10, 0.95, 1.23, 1.13, 0.93, 1.67, 2.51, 0.98, 0.80, 0.985282715, 1.639191663, 1.382705136, 0.48874792, 1.398489625, 1.265029243, 1.286735024, 1.146877561, 1.827147032, 2.534746454, 1.882753246, 1.475607608, 1.863883108, 2.277982827, 3.540565615, 1.707320794, 2.568958439, 14.54490887, 3.237073047, 3.7630828, 3.942674884, 1.578639917, 2.996229014, 1.089460882, 1.154768444, 1.979792303, 1.357551379)
  
  df <- right_join(data, pred) %>% 
    mutate(upper = index + U,
           lower = index - L) %>%
    tbl_df()
  
  ggplot(df, aes(year, index)) + 
    geom_point(colour = "grey") +
    geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey") + 
    labs(x = "Year", y = "Egg Deposition (trillions)\n") +
    if(fit) geom_line(aes(year, pred), alpha = 0.8) 
  
}


# ---------------------------------------------------------------------------- #
# PLOTS FOR DATA SECTION
# ---------------------------------------------------------------------------- #
# prefix with d for data
# d1 » catch time series
# d2 » spawn weight-at-age
# d3 » commercial weight-at-age
# d4 » commercial age-proportions
# d5 » spawn sample age-proportions
# d6 » Survey- and model-estimated egg deposition
d1 <- plot.catch(D, nm = "data_ct_raw", y = "Catch (tons)")
d2 <- plot.waa(D, nm = "data_sp_waa")
d3 <- plot.waa(D, nm = "data_cm_waa")
d4 <- plot.comp(D, nm = "data_cm_comp")
d5 <- plot.comp(D, nm = "data_sp_comp")
d6 <- plot.eggdepfit(D, sfx = "egg_dep", fit = TRUE)
plot.ssb(D)

plot.ssb.post(D)


# Run models:

# clean up directory - remove unnecessary files to run
# add ".pin" when wanting to save .pin file, or any other file to save
setwd("HER_2018forec/")
need_files <- c(".tpl", ".ctl", ".dat", ".R", ".csv")
files_present <- list.files()
keep_index <- unlist(sapply(1:length(need_files), function(x) which(grepl(need_files[x], files_present))))
rm_files <- files_present[-keep_index]
remove <- sapply(1:length(rm_files), function(x) unlink(rm_files[x], TRUE))
# setwd("..") # backs out one folder

# LS model 2018
LS_year <- read_csv("LS_2018forec_results_byyear.csv")

## compile model
setup_admb()
compile_admb("her", verbose=TRUE)

## run MAP
run_admb("her")

years <- 1980:2017
ages <- 3:8

#Maturity
readMat("mat", file="her.rep", nrow = length(years))
readVec("mat_params[1]", file="her.rep")
readVec("mat_params[2]", file="her.rep")

# Natural mortality
natmat <- readMat("Mij", file = "her.rep", nrow = length(years))[,1]
plot(natmat ~ years, ylim=c(0, max(natmat)*1.1), type="l", lwd=2)

# Selectivity
sel <- readMat("Sij", file = "her.rep", nrow = length(years))[,1]
plot(sel ~ years, ylim=c(0, max(natmat)*1.1), type = "l", lwd = 2)

## read report from initial MAP run - maximum a posteriori estimation (i.e.
## maximum likelihood using priors!)
ssb <- readVec("ssb", file="her.rep")
## put in thousands 
# ssb <- ssb/1000

plot(ssb ~ years, ylim = c(0, max(ssb)*1.1), type = "l", lwd = 2)
lines(LS_year$tot_sp_B_tons ~ years, type = "l", col = "red", lwd = 2)

## run simulation with seed 123
run_admb("her", extra.args="-sim 123")

## run MCMC
run_admb("her", extra.args="-mcmc 1000 -mcsave 10")
run_admb("her", extra.args="-mceval")

## posterior distributions
ssb_ps <- read.table("ssb.ps")
natural_ps <- read.table("natural.ps", header=TRUE)

par(mfrow=c(2,1))
plot(natural_ps[,1])
abline(h=median(natural_ps[,1]), col="red")
plot(natural_ps[,2])
abline(h=median(natural_ps[,2]), col="red")

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