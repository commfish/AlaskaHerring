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

# Functions for getting BRPs  ----

# Derive steepness (h) from compensation ratio (reck)
get_steepness <- function(rec_mod) {
  
  # Ricker
  if(rec_mod == 1) { h <- reck / (4 + reck) }
  
  # Beverton-Holt
  if(rec_mod == 2) { h <- 0.2 * reck ^ (4 / 5) }
  
  return(h)
}

# Survivorship
get_survivorship <- function(fmort) {
  
  lx <- matrix(ncol = length(unique(age)), 
               nrow = length(unique(fmort)))
  
  a_plus <- max(length(age)) 
  
  # Loop over all input values of mortality
  for(f in 1:length(fmort)) {
    
    # Initialize matrix
    lx[f,1] <- 1  
    
    for(a in 2:(a_plus-1)) {
      lx[f,a] <- lx[f,a-1] * exp(-mort[a-1] - fmort[f] * sel[a-1]) }
    
    # Plus group
    lx[f,a_plus] <- (lx[f,a_plus-1] * exp(- mort[a_plus-1] - fmort[f] * sel[a_plus-1])) / (1 - exp(- mort[a_plus] - fmort[f] * sel[a_plus]))    
  }
  return(lx)
}

# Yield-per-recruit (note that ypr = phi_q * fmort in Martell et al. 2009)
get_ypr <- function(fmort, lx) { 
  
  ypr <- matrix(ncol = length(unique(age)), 
                nrow = length(unique(fmort))) 
  
  # Loop over all input values of mortality
  for(f in 1:length(fmort)) {
    for(a in 1:length(age)) {
      ypr[f,a] <- wa[a] * (fmort[f] * sel[a]) / (fmort[f] * sel[a] + mort[a]) * 
        lx[f,a] * (1 - exp(-(fmort[f] * sel[a] + mort[a]))) 
    }
  }
  
  # Sum across ages for each F value
  ypr <- rowSums(ypr) 
  return(ypr)
}

# Vulnerable biomass-per-recruit (phi_b in Martell et al. 2009)
get_bpr <- function(fmort, lx) {
  
  bpr <- matrix(ncol = length(unique(age)), 
                nrow = length(unique(fmort)))
  
  for(f in 1:length(fmort)) {
    for(a in 1:length(age)) {
      bpr[f,a] <- lx[f,a] * wa[a] * sel[a]
    }
  }
  # Sum across ages for each F value
  bpr <- rowSums(bpr)
  return(bpr)
}

# Function to get spawner-per-recruit (phi_e in Martell et al. 2009)
get_spr <- function(fmort, lx) { 
    
    # spr matrix of same dimensions as N matrix
    spr <- matrix(ncol = length(unique(age)), 
                  nrow = length(unique(fmort)))
    
    # Loop over all input values of fishing mortality
    for(f in 1:length(fmort)) {
      for(a in 1:length(unique(age))){
        spr[f,a] <- wa[a] * mat[a] * lx[f,a] }}
    
    # Sum across ages for each F value
    spr <- rowSums(spr) 
    return(spr)
  }

# Derive stock-recruitment parameters 
get_recruit_pars <- function() {
  
  # Ricker model
  if(rec_mod == 1) {
    alpha <- exp(log(5*h)/0.8)/spr0
    beta <- log(5*h)/(0.8*spr0*r0) }
  
  # Beverton-Holt model
  if(rec_mod == 2) {
    alpha <- spr0*(1-h)/(4*h)
    beta <- (5*h-1)/(4*h*r0) }
  
  recruit_params <- list(alpha, beta)
  return(recruit_params)
}

# Get recruitment conditioned on alpha/beta and spr
get_recruitment <- function(spr) { 
    # Ricker model (Eqn 10a, put in terms of spr)
  if(rec_mod == 1) {
    r <- log(alpha*spr)/(beta*spr) }
  
  # Beverton-Holt model (Eqn 9)
  if(rec_mod == 2) {
    r <- (spr-alpha)/(beta*spr) }
  
  return(r)
}

# Get equilibrium conditions
get_equilibrium <- function(ff) {
  
  lx_e <- get_survivorship(fmort = ff)  
  ypr_e <- get_ypr(fmort = ff, lx = lx_e)
  spr_e <- get_spr(fmort = ff, lx = lx_e) 
  bpr_e <- get_bpr(fmort = ff, lx = lx_e)
  r_e <- get_recruitment(spr = spr_e)
  
  out <- NULL
  out$ypr <- ypr_e
  out$spr <- spr_e
  out$r <- r_e
  out$bpr <- bpr_e
  out$b <- bpr_e * r_e
  out$s <- spr_e * r_e
  out$s_ratio <- (spr_e * r_e) / s0 # Express spawners as a proportion of unfished spawners
  out$y <- ypr_e * r_e

  return(out)
}

# Derivative of ypr with respect to F (use uniroot() to solve for dY(F)/df=0)
dfx.dx <- function(ff, delta) {
  y1 <- get_equilibrium(ff = ff-delta/2)$y
  y2 <- get_equilibrium(ff = ff+delta/2)$y
  approx.gradient <- (y2-y1)/delta
  return(approx.gradient)
}

# Get F crash (where F is minimized and S ~ 0 using the bisection method - MK helped me with this.
bisect <- function(Fmin = 1, Fmax = 100){
  for(b in 1:100000){
    F_tst <- (Fmin + Fmax)/2 # update
    s_tmp <- get_equilibrium(ff = F_tst)$s
    if(round(s_tmp,0) == 0 & (Fmax - Fmin) > 0.0001){ return(F_tst)
    } else if(round(s_tmp,4) > 0) { Fmin <- F_tst 
    } else if(round(s_tmp,4) < 0) { Fmax <- F_tst }
  }
  print('max iter')
}

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

