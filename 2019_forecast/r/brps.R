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
TIME_PERIOD <- paste0("1980-",YEAR-1)

# Number of interations (MCMC samples) and thinning rate (if 10, saves every
# 10th sample) - make sure these values are the same as in her.r
niter <- 11110000 # this value give you nout=10000
thin <- 1000
burn_in <- 0.1 # 10% (this is the first 10% of the thinned/saved iterations)
nout <- round((niter / thin + 1) - burn_in * (niter / thin + 1), 0)

source(paste0(YEAR, "_forecast/r/helper.r"))

# Directories ----
main_dir <- getwd()

# Directory for biological reference point results for a given model and time
# period of data
brp_dir <- file.path(main_dir, paste0(YEAR, "_forecast/results/reference_points"))
dir.create(brp_dir, showWarnings = FALSE)
brp_dir <- file.path(brp_dir, paste0(MODEL_VERSION))
dir.create(sub_dir, showWarnings = FALSE)
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
                     burn = burn_in * (niter / thin + 1)) {
  
  require(data.table)
  df <- fread(paste0(fn, ".ps"))
  colnames(df) <- paste(sage:nage)
  df[, Year := rep(syr:lyr, n)]
  df[, iter := rowid(Year)] 
  df <- melt(df, id.vars = c("Year", "iter"), variable.name = "Age")
  df <- df[iter > burn, ] # Eliminate burn in
  
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

# Spawner weight-at-age
wa <-  as.data.frame(D[["data_sp_waa"]])
colnames(wa) <- c("Year", paste(sage:nage))
wa <-	wa %>% 
  melt(id.vars = c("Year"), variable.name = "Age", value.name = "weight") %>% 
  filter(Year >= syr) %>% 
  group_by(Age) %>% 
  dplyr::summarise(weight = mean(weight)) %>% 
  pull(weight)

# Posterior samples for r0 and reck
pars <- data.table(D[["post.samp"]][,4:5]) # (see her.ctl to check parameter order)
colnames(pars) <- c("log_r0", "log_reck")
pars[, iter := .I] # row number = iteration number
burn <- round(burn_in * (niter / thin + 1), 0) # index for end of burn in
pars <- pars[burn:(niter/thin),] # Remove burn-in

# Functions for getting BRPs

# Inputs
age <- sage:nage
fmort <- seq(0, 1, 0.01)
rec_mod <- 1 # 1 = Ricker, 2 = Beverton-Holt
m_iter = 1112

# Natural mortality
mort <- survival[iter == m_iter, ]
mort[, mort := -log(survival)]
mort <- mort$mort

# Selectivity
sel <- selectivity[iter == m_iter, ]
sel <- sel$selectivity

# Maturity
mat <- maturity[iter == m_iter, ]
mat <- mat$maturity

# Recruitment k (compensation ratio) and Equilibrium recruitment
pars_s <- pars[iter == m_iter, ]
reck <- exp(pars_s$log_reck) 
r0 <- exp(pars_s$log_r0)  

# lx = survivorship under unfished conditions  
get_survivorship <- function() {
  
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

# Yield-per-recruit
get_ypr <- function() { 
  
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

# Function to get spawner-per-recruit
get_spr <-
  function() { 
    
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

# Derive steepness (h) from compensation ratio (reck)
get_steepness <- function() {
  
  # Ricker
  if(rec_mod == 1) { h <- reck / (4 + reck) }
  
  # Beverton-Holt
  if(rec_mod == 2) { h <- 0.2 * reck ^ (4 / 5) }
  
  return(h)
}

# Derive stock-recruitment parameters 
get_recruit_pars <- function() {
    
    # Beverton-Holt model
    if(rec_mod == 1) {
      alpha <- spr0*(1-h)/(4*h)
      beta <- (5*h-1)/(4*h*r0) }
    
    # Ricker model
    if(rec_mod == 2) {
      alpha <- exp(log(5*h)/0.8)/spr0
      beta <- log(5*h)/(0.8*spr0*r0) }
    
    recruit_params <- list(alpha, beta)
    return(recruit_params)
  }

# Get recruitment conditioned on alpha/beta and spr
get_recruitment <- function() { 
  # Beverton-Holt model (Eqn 9)
  if(rec_mod == 1) {
    r <- (spr-alpha)/(beta*spr) }
  
  # Ricker model (Eqn 10a, put in terms of spr)
  if(rec_mod == 2) {
    r <- log(alpha*spr)/(beta*spr) }
  
  return(r)
}

h <- get_steepness()

get_equilibrium <- function(fmort = fmort) {
  
  lx <- get_survivorship()  
  ypr <- get_ypr()
  spr <- get_spr() 
  spr0 <- spr[1] # unfished spawner per recruit
  alpha <- get_recruit_pars()[[1]]
  beta <- get_recruit_pars()[[2]]
  r <- get_recruitment()
  s0 <- r0 * spr0 # unfished spawning biomass 
  
  out <- NULL
  out$ypr <- ypr
  out$spr <- spr
  out$spr0 <- spr0
  out$h <- h
  out$alpha <- alpha
  out$beta <- beta
  out$r <- r
  out$s0 <- s0
  out$s <- spr*r
  out$s_ratio <- (spr*r)/s0 # Express spawners as a proportion of unfished spawners
  out$y <- ypr*r
  
  return(out)
}

equilibrium <- get_equilibrium(fmort = fmort)


# Function to get derivative of ypr with respect to F (use uniroot() to solve for dY(F)/df=0)
dfx.dx <- function(fmort, get_equilibrium, delta, h) {
  y1 <- get_equilibrium(fmort = fmort-delta/2)$y
  y2 <- get_equilibrium(fmort = fmort+delta/2)$y
  approx.gradient <- (y2-y1)/delta
  return(approx.gradient)
}

# Get F crash (where F is minimized and S ~ 0 using the bisection method - MK helped me with this.
bisect <- function(Fmin = 0, Fmax = 1){
  for(b in 1:10000){
    F_tst <- (Fmin + Fmax)/2 # update
    s_tmp <- get_equilibrium(fmort = F_tst)$s
    if(round(s_tmp,4) == 0 & (Fmax - Fmin) > 0.0002){ return(F_tst)
    } else if(round(s_tmp,4) > 0) { Fmin <- F_tst 
    } else if(round(s_tmp,4) < 0) { Fmax <- F_tst }
  }
  print('max iter')
}

mod_names <- c("Beverton-Holt", "Ricker")

# Get biological reference points 
brps <- data.frame(Model = NA, Fmsy = NA, MSY = NA, Fcrash = NA)

# Loop over the 3 S-R models
for(rec_mod in 1:2){
  # uniroot will only take one value at a time, it will pick values on the interval and feed them as FF's to dfx.dx
  Fmsy <- uniroot(f = dfx.dx, interval = c(0.01,1), tol = 0.000001, 
                  get_equilibrium = get_equilibrium, delta = 0.001)$root[1]
  msy <- get_equilibrium(fmort = Fmsy)$y
  Fcrash <- bisect()
  brps[mod,] <- c(mod_names[mod], round(Fmsy, 4), round(msy, 4), round(Fcrash, 4))
}


for(MOD in 1:3){
  
  out <- yield_f(FF, h, MOD = MOD)
  refs <- brps[MOD,]
  
  par(mar=c(4,4,3,3), mfrow = c(2,2))
  
  plot(out$s_ratio, out$r, type = "l", xlab = "S(F)", ylab = "R(F)",
       xlim = c(0, max(out$s_ratio)), ylim = c(0, max(out$r)*1.1))
  abline(0,1, col = "grey", lty = 2)
  
  plot(out$s_ratio, out$y, type = "l", xlab = "S(F)", ylab = "Y(F)",
       xlim = c(0, max(out$s_ratio)), ylim = c(0, max(out$y)*1.1))
  
  text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
       line2user(line=2, side=3), MOD_NAMES[MOD], xpd=NA, cex=2, font=2)
  
  plot(FF, out$y, type = "l", xlab = "F", ylab = "Y(F)",
       xlim = c(0, as.numeric(refs$Fcrash)*1.2), ylim = c(0, max(out$y)*1.1))
  abline(h = refs$MSY, col = "grey", lty = 2)
  abline(v = refs$Fmsy, col = "grey", lty = 2)
  abline(v = refs$Fcrash, col = "grey", lty = 2)
  
  plot(FF, out$s_ratio, type = "l", xlab = "F", ylab = "S(F)",
       xlim = c(0, as.numeric(refs$Fcrash)*1.2), ylim = c(0, max(out$s_ratio)*1.1))
  abline(v = refs$Fcrash, col = "grey", lty = 2)
  
}
