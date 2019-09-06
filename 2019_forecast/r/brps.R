# Biological reference points
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov
# Last edited: 2019-08-11

# Estimates equilibrium biological reference points using posterior samples of
# R0 (equilibrium recruitment), maturity, natural mortality, and reck
# (recruitment compensation ratio).

# User-defined variables ----
YEAR <- 2019 # Forecast year
MODEL_VERSION <- "HER_best_condEffort.12_322"   # HER with best HER parameterization by AIC, conditioned on effort
# MODEL_VERSION <- "HER_best_condCatch.12_322"   # HER with best HER parameterization by AIC, conditioned on effort
TIME_PERIOD <- paste0("1980-",YEAR-1)

# Number of interations (MCMC samples) and thinning rate (if 10, saves every
# 10th sample) - make sure these values are the same as in her.r
niter <- 11110000 # this value give you nout=10000
thin <- 1000
burn_in <- 0.1 # 10% (this is the first 10% of the thinned/saved iterations)
nout <- round((niter / thin + 1) - burn_in * (niter / thin + 1), 0)

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
source(paste0(YEAR, "_forecast/r/helper.r"))
# install.packages("svMisc")
library(svMisc)

# Directories ----
main_dir <- getwd()

# Directory for biological reference point results for a given model and time
# period of data
brp_dir <- file.path(main_dir, paste0(YEAR, "_forecast/results/reference_points"))
dir.create(brp_dir, showWarnings = FALSE)
brp_dir <- file.path(brp_dir, paste0(MODEL_VERSION))
dir.create(brp_dir, showWarnings = FALSE)
brp_dir <- file.path(brp_dir, paste0(TIME_PERIOD))
dir.create(brp_dir, showWarnings = FALSE)

# Read in data and model posterior samples ----
setwd(file.path(main_dir, paste0(YEAR, "_forecast/admb/", MODEL_VERSION)))

D <- read_admb("her")

# Model dimensions
sage <- D[['sage']]
nage <- D[['nage']]
ages <- sage:nage

syr <- D[['mod_syr']]
nyr <- D[['mod_nyr']]
yrs <- syr:nyr


# Function to read time-varying parameter posterior samples
read_tvpar_ps <- function(fn = "maturity", 
                     syr = D[["mod_syr"]], 
                     lyr = D[["mod_nyr"]], 
                     n = niter / thin + 1,
                     burn = burn_in * (niter / thin + 1)
                     # burn = burn_in * (niter / thin)
                     ) {
  
  require(data.table)
  df <- fread(paste0(fn, ".ps"))
  colnames(df) <- paste(sage:nage)
  df[, Year := rep(syr:lyr, n)]
  df[, iter := rowid(Year)] 
  df <- melt(df, id.vars = c("Year", "iter"), variable.name = "Age")
  df <- df[iter >= burn, ] # Eliminate burn in
  # df <- df[iter < max(iter)]
  
  # Selectivity has a row for the forecast year. Get rid of that row before
  # getting mean for each iter
  if(fn == "selectivity") {
    df <- df[Year != YEAR]
  }
  # akin to dplyr::summarize ~ get mean value for each iteration (mean over all
  # years)
  df <- df[, list(value = mean(value)),
           by = .(iter, Age)]
  colnames(df) <- c("iter", "Age", paste0(fn))
  df[, age := as.integer(as.character(Age))]
  return(df)
}

# Maturity, selectivity, survival (exp(-M))
maturity <- read_tvpar_ps(fn = "maturity")
selectivity <- read_tvpar_ps(fn = "selectivity", lyr = D[["mod_nyr"]] + 1)
survival <- read_tvpar_ps(fn = "survival")

# Posterior samples for r0 and reck
pars <- data.table(D[["post.samp"]][,4:5]) # (see her.ctl to check parameter order)
colnames(pars) <- c("log_r0", "log_reck")
pars[, iter := .I] # row number = iteration number
burn <- round(burn_in * (niter / thin + 1), 0)  # index for end of burn in
pars <- pars[burn:(niter/thin), ] # Remove burn-in

# Spawner weight-at-age
wa <-  as.data.frame(D[["data_sp_waa"]])
colnames(wa) <- c("Year", paste(sage:nage))
wa <-	wa %>% 
  melt(id.vars = c("Year"), variable.name = "Age", value.name = "weight") %>% 
  filter(Year >= syr) %>% 
  group_by(Age) %>% 
  dplyr::summarise(weight = mean(weight)) %>% 
  pull(weight)

# Run analyis ----
# Model names (currently only uses Ricker)
mod_names <- c("Ricker", "Beverton-Holt")

age <- sage:nage
fmort <- seq(0, 20, 0.01)
rec_mod <- 1 # 1 = Ricker, 2 = Beverton-Holt

output <- list()
brps <- list()

for(i in 1:length(unique(pars$iter))) {
  
  # Inputs
  # Natural mortality
  mort <- survival[iter == unique(iter)[i], ]
  mort[, mort := -log(survival)]
  mort <- mort$mort
  
  # Selectivity
  sel <- selectivity[iter == unique(iter)[i], ]
  sel <- sel$selectivity
  
  # Maturity
  mat <- maturity[iter == unique(iter)[i], ]
  mat <- mat$maturity
  
  # Recruitment k (compensation ratio) and Equilibrium recruitment
  pars_s <- pars[iter == unique(iter)[i], ]
  reck <- exp(pars_s$log_reck) + 1 # FLAG - in her.tpl, he adds 1. Assuming we should do the same.
  r0 <- exp(pars_s$log_r0)  
  h <- get_steepness(rec_mod = rec_mod)
  
  # Get BRPs

  # Unfished survival, spawning biomass per recruit, and unfished spawning biomass
  lx0 <- get_survivorship(fmort = 0)
  spr0 <- get_spr(fmort = 0, lx = lx0)
  bpr0 <- get_bpr(fmort = 0, lx = lx0)
  s0 <- r0 * spr0
  
  # Recruitment parameters
  alpha <- get_recruit_pars()[[1]]
  beta <- get_recruit_pars()[[2]]
  
  # uniroot will only take one value at a time, it will pick values on the interval and feed them as FF's to dfx.dx
  Fmsy <- uniroot(f = dfx.dx, interval = c(0.01,100), tol = 0.000001, 
                  delta = 0.0001)$root[1]
  msy <- get_equilibrium(ff = Fmsy)$y
  # Fcrash <- bisect()
  
  # Get components of B0
  lx_Fmsy <- get_survivorship(fmort = Fmsy)
  ypr_Fmsy <- get_ypr(fmort = Fmsy, lx = lx_Fmsy) # Fmsy * phi_q evaluated at Fmsy
  spr_Fmsy <- get_spr(fmort = Fmsy, lx = lx_Fmsy) # phi_e evaluated at Fmsy
  r_Fmsy <- get_recruitment(spr = spr_Fmsy)
  
  bmsy <- r_Fmsy * spr_Fmsy  # Equilibrium spawning biomass at Fmsy
  b0 <- r0 * spr0 # Equilibrium unfished spawning biomass
  
  brps[[i]] <- data.frame(Model = mod_names[rec_mod], iter = i, h = round(h, 4), 
                     alpha = round(alpha, 3), beta = round(beta, 5), s0 = round(s0, 0), 
                     spr0 = round(spr0, 0), Fmsy = round(Fmsy, 4), MSY = round(msy, 0), 
                     #Fcrash = round(Fcrash, 4), 
                     Bmsy = round(bmsy, 0), B0 = round(b0, 0))
  
  out <- get_equilibrium(fmort)
  
  out <- out %>% 
    map(~data.frame(.x)) %>% 
    imap(~mutate(.x, variable = .y)) %>%
    map(~mutate(.x, iter = i))
  out <- bind_rows(out)
  colnames(out) <- c("value", "variable", "iter")
  output[[i]] <- out
  
  # progress(value = i, progress.bar = TRUE, init = 1)
  # Sys.sleep(0.03)
  # if (i == max(length(unique(pars$iter)))) cat("Done!\n")
}

# Summarize/save output ----

brps <- bind_rows(brps)
output <- bind_rows(output)

brps <- melt(brps, id.vars = c("Model", "iter"))

brps <- data.table(brps)
brps[, `:=` (mean = mean(value),
             median = median(value),
             # 95% 
             q025 = quantile(value, 0.025),
             q975 = quantile(value, 0.975),
             # 50% 
             q250 = quantile(value, 0.250),
             q750 = quantile(value, 0.750)),
     by = .(variable)] 

brp_sum <- unique(brps, by = c("variable", "mean", "median", "q025", "q975", "q250", "q750"))
  
# hist(brps[variable == "B0"]$value)

output <- data.table(output)
output[, fmort := rep(fmort, length(unique(iter))*length(unique(variable)))]
# akin to dplyr::summarize()
output <- output[, list(mean = mean(value),
                        median = median(value),
                        # 95% 
                        q025 = quantile(value, 0.025),
                        q975 = quantile(value, 0.975),
                        # 50% 
                        q250 = quantile(value, 0.250),
                        q750 = quantile(value, 0.750)),
                 by = .(variable, fmort)]

write_csv(brps, paste0(brp_dir, "/raw_brps.csv"))
write_csv(brp_sum, paste0(brp_dir, "/brps_sum.csv"))
write_csv(output, paste0(brp_dir, "/raw_output.csv"))

# ---
brps <- read_csv(paste0(brp_dir, "/raw_brps.csv"))
brp_sum <- read_csv(paste0(brp_dir, "/brps_sum.csv"))
output <- read_csv(paste0(brp_dir, "/raw_output.csv"))
brps <- data.table(brps)
brp_sum <- data.table(brp_sum)
output <- data.table(output)

# MLE estimates ----
names(D)


mat <- read_D("mat")
# mat <- LS_byage %>%
#   select(Age, maturity) %>%
#   group_by(Age) %>%
#   dplyr::summarize(out = mean(maturity)) %>%
#   pull(out)
# 
# sel <- LS_byage %>%
#   select(Age, maturity) %>%
#   group_by(Age) %>%
#   dplyr::summarize(out = mean(maturity)) %>%
#   pull(out)
  
sel <- read_D("Sij")
mort <- read_D("Mij")
r0 <- exp(D[["theta"]][4])
reck <- exp(D[["theta"]][5]) + 1 # FLAG - in her.tpl, he adds 1. Assuming we should do the same.

h <- get_steepness(rec_mod = rec_mod)

# Unfished survival, spawning biomass per recruit, and unfished spawning biomass
lx0 <- get_survivorship(fmort = 0)
spr0 <- get_spr(fmort = 0, lx = lx0)
bpr0 <- get_bpr(fmort = 0, lx = lx0)
s0 <- r0 * spr0

# Get BRPs
out <- get_equilibrium(ff = fmort)

# Recruitment parameters
alpha <- get_recruit_pars()[[1]]
beta <- get_recruit_pars()[[2]]

# uniroot will only take one value at a time, it will pick values on the interval and feed them as FF's to dfx.dx
Fmsy <- uniroot(f = dfx.dx, interval = c(0.01,100), tol = 0.000001, 
                delta = 0.0001)$root[1]
msy <- get_equilibrium(ff = Fmsy)$y
# Fcrash <- bisect()

# Get components of B0
lx_Fmsy <- get_survivorship(fmort = Fmsy)
ypr_Fmsy <- get_ypr(fmort = Fmsy, lx = lx_Fmsy) # Fmsy * phi_q evaluated at Fmsy
spr_Fmsy <- get_spr(fmort = Fmsy, lx = lx_Fmsy) # phi_e evaluated at Fmsy
r_Fmsy <- get_recruitment(spr = spr_Fmsy)

bmsy <- r_Fmsy * spr_Fmsy  # Equilibrium spawning biomass at Fmsy
b0 <- r0 * spr0

notransformation <- b0*.907
transformreck <- b0*.907
.25*b0*.907

spb_sum <- ps_byyear(save = FALSE, fn = "sp_B", unit_conv = 0.90718)
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
  geom_hline(yintercept = transformreck, lty = 2, col = "red") +
  geom_hline(yintercept = notransformation, lty = 2, col = "blue") +
  scale_shape_manual(values = 1) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  scale_y_continuous(limits = c(0, max(LS_byyear$surv_est_spB)), labels = scales::comma) +
  labs(x = NULL, y = "Spawning biomass (tons)\n", shape = NULL) +
  theme(legend.position = c(0.25, 0.8))

# Graphics ----
y <- output[variable == "y", ]

ggplot(y) +
  geom_line(aes(x = fmort, y = mean)) +
  geom_line(aes(x = fmort, y = median), col = "red") +
  geom_ribbon(aes(x = fmort, ymin = q025, ymax = q975),
              alpha = 0.4, fill = "grey80") +
  geom_ribbon(aes(x = fmort, ymin = q250, ymax = q750),
              alpha = 0.4, fill = "grey60") +
  geom_hline(yintercept = c(brp_sum[variable == "MSY",]$median)) +
  coord_cartesian(ylim = c(0, max(y$q975)), xlim = c(0, 7)) +
  labs(x = "Fishing mortality", y = "Yield")
  
hist(brps[variable == "B0"]$value)
brps[variable == "Bmsy"]
brps[variable == "Fmsy"]
(15234 / 0.90718) / 2


ggplot(output, aes)


# 


refs <- brps[rec_mod,]
refs <- brps[[i]]

par(mar=c(4,4,3,3), mfrow = c(2,2))
plot(fmort, out$s, type = "l")
abline(h = refs$Bmsy, col = "grey", lty = 2)
abline(v = refs$Fmsy, col = "grey", lty = 2)


# abline(h = refs$B0, col = "grey", lty = 2)
# abline(v = 0, col = "grey", lty = 2)
# abline(h = 0, col = "grey", lty = 2)
# abline(v = refs$Fcrash, col = "grey", lty = 2)
#
#
plot(fmort, out$b, type = "l")
# abline(h = refs$Bmsy, col = "grey", lty = 2)
# abline(v = refs$Fmsy, col = "grey", lty = 2)
# abline(h = refs$B0, col = "grey", lty = 2)
#
plot(out$s_ratio, out$r / r0, type = "l", xlab = "S(F)", ylab = "R(F)",
     xlim = c(0, max(out$s_ratio)), ylim = c(0, max(out$r / r0)*1.1))
# abline(0,1, col = "grey", lty = 2)
#
plot(out$s_ratio, out$y, type = "l", xlab = "S(F)", ylab = "Y(F)",
     xlim = c(0, max(out$s_ratio)), ylim = c(0, max(out$y)*1.1))
#
# text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
#      line2user(line=2, side=3), mod_names[rec_mod], xpd=NA, cex=2, font=2)
#
plot(fmort, out$y, type = "l", xlab = "F", ylab = "Y(F)",
     xlim = c(0, 30), ylim = c(0, max(out$y)*1.1))
# abline(h = refs$MSY, col = "grey", lty = 2)
# abline(v = refs$Fmsy, col = "grey", lty = 2)
# abline(v = refs$Fcrash, col = "grey", lty = 2)
#
# plot(fmort, out$s_ratio, type = "l", xlab = "F", ylab = "S(F)",
#      xlim = c(0, as.numeric(refs$Fcrash)*1.2), ylim = c(0, max(out$s_ratio)*1.1))
# abline(v = refs$Fcrash, col = "grey", lty = 2)
# abline(v = refs$Fmsy, col = "grey", lty = 2)

#}

# Need more details from SM before I'd use these.
# # equation (8) in Martell et al. 2009
# if(rec_mod == 1) {
#   b0 <- - (log(reck) * bpr0 * spr_Fmsy * msy) / (log(spr0/(reck*spr_Fmsy)) * spr0 * ypr_Fmsy * Fmsy)
# }
# # equation (5) in Martell et al. 2009
# if(rec_mod == 2) {
#   b0 <- (msy * bpr0 * (reck - 1)) / (ypr_Fmsy * (reck - spr0/spr_Fmsy))
# }

P <- read_fit("her")

obj_fn <- P[["nlogl"]]

data.frame(like = c("sp_comp", "cm_comp", "egg_dep", 
                    "milt", "sr", "catch",
                    "rinit", "rbar", "fbar", "mdevs",
                    "prior_logM", "prior_logrinit", "prior_logrbar",
                    "prior_logro", "prior_logreck", "prior_precr"), 
           nll = c(D[["nll"]], D[["penll"]], D[["calcPriors()"]])) %>% 
  # mutate(ll = -1*nll,
  #        like = round(exp(-1*nll),0))
  mutate(perc = nll/obj_fn) %>% 
  summarise(sum(perc))

sum(c(D[["nll"]], D[["penll"]], D[["calcPriors()"]]))
sum(D)
names(D)
D[["calcPriors()"]]
