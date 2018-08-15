
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

# in SPAWN folder, numbers from excel spread sheet for spawn deposition.
LS_byyear$surv_est_spB <- c(35000,30000,29500,23500,38500,31000,25000,46000,58500,27000,23000,23500,
                             43351,37150,14941,34990,40827,28611,34942,44554,57988,58756,40366,55769,
                             69907,101305,66111,84501,247088,110946,126230,161904,62518,103267,48561, 
                             58183,77973,46919)

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

#Print out natural mortality
df <- data.frame(D[["year"]], 
                 D[["Mij"]])
colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
write.table(df,sep="\t", row.names=FALSE, col.names=FALSE)

#Print out maturity
df <- data.frame(D[["year"]], 
                 D[["mat"]])
colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
write.table(df,sep="\t", row.names=FALSE, col.names=FALSE)

#Print out selectivity
# number of yrs in the model 
nyr <- D[["mod_nyr"]] - D[["mod_syr"]] + 1

df <- data.frame(D[["year"]], 
                 D[["Sij"]][1:nyr, ]) # change number of rows
colnames(df) <- c("Year",paste(D[['sage']]:D[['nage']]))
write.table(df,sep="\t", row.names=FALSE, col.names=FALSE)

# 1. Time series of mature biomass (pre-fishery) and spawning biomass (post-fishery)
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
  geom_point(data = df, aes(x = Year, y = tons, colour = Model, linetype = Model, shape = Model)) +
  geom_line(data = df, aes(x = Year, y = tons, colour = Model, linetype = Model, shape = Model)) +
  geom_point(data = srv_index, aes(x = Year, y = srv_spB), shape = "*", size = 5) +
  #scale_colour_grey() +
  theme(legend.position = c(0.1, 0.8)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  labs(x = "", y = "Spawning biomass (short tons)\n")

# 2. Recruitment
df <- data.frame(Year = seq(D[["mod_syr"]] + D[["sage"]],
                      D[["mod_nyr"]] + 1, 1),
           recruits = D[["recruits"]],
           Model = "HER") %>% 
  bind_rows(LS_byyear %>% 
  select(Year = year, recruits = init_age_3) %>% #
  mutate(Model = "LS"))

ggplot(df, aes(x = Year, y = recruits, colour = Model, shape = Model)) +
  geom_point() +
  geom_line() +
  # scale_colour_grey() +
  theme(legend.position = c(0.1, 0.8)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
  labs(x = "", y = "Age-3 recruits (millions)\n")

# 3. Egg deposition
df <- as.data.frame(D[["data_egg_dep"]][10:47, ])#[seq(tot_yrs - nyr + 1, tot_yrs, 1), ]) 
colnames(df) <- c("Year", "obs", "log_se")
df %>% filter(Year >= D[["mod_syr"]]) %>% 
  bind_cols(data.frame(HER = D[["pred_egg_dep"]],
                       LS = LS_byyear$tot_est_egg[1:38]))%>% 
  gather("Model", "trillions", -c(Year, obs, log_se)) -> df

ggplot(df, aes(x = Year)) +
  geom_point(aes(y = obs), shape = "*", size = 5) +
  geom_line(aes(y = trillions, colour = Model,
                linetype = Model, shape = Model)) +
  geom_point(aes(y = trillions, colour = Model,
                linetype = Model, shape = Model)) +
  # scale_colour_grey() +
  theme(legend.position = c(0.1, 0.8)) +
  scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +

  labs(x = "", y = "Eggs spawned (trillions)\n")

# 4. spawn age comps, cast net
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
  mutate(Predicted = "HER",
         Age = factor(Age, levels = c("3", "4", "5", "6", "7", "8"),
                      labels = c("3", "4", "5", "6", "7", "8+"))) %>% 
  bind_rows(LS_byage %>% 
              select(Year, Age, proportion = spawnage_comp_est) %>% 
              mutate(Predicted = "LS") ) -> pred

ggplot() +
  geom_bar(data = obs, aes(x = Age, y = proportion), 
           stat = "identity", colour = "lightgrey", fill = "lightgrey") +
  geom_line(data = pred, aes(x = Age, y = proportion, colour = Predicted, group = Predicted), size = 1) +
  facet_wrap(~ Year, dir = "v", ncol = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)) +
  labs(x = '\nAge', y = 'Proportion-at-age\n') 



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

## from tech doc:
## 1) first fit to sitka data
## 2) save .par file as her.pin
D$