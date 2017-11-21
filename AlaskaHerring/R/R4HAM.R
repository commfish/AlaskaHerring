# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# R4HAM.R 
# Author: Steven Martell
# Sep 27, 2016
# Objective:
# 	This is a simple example of how I structure my R-scripts so they are
#   easy to maintain and portable.
#
#
#
# ---------------------------------------------------------------------------- #

#
# LIBRARIES
#
library(ggplot2)
library(dplyr)
library(tidyr)


#
#	OTHER R-SRCIPTS
#
source(file.path("globals.R"))
source(file.path("plotIndex.R"))
source(file.path("ham.R"))


# ---------------------------------------------------------------------------- #
# READ MODEL OUTPUTS AND STORE IN A LIST
# ---------------------------------------------------------------------------- #
# MODELDIRS will recursively look in all folders specified in the file path.
"~/Library/Mobile Documents/com~apple~CloudDocs/SSI/ADFG_HerringModel/R"
.MODELDIRS <- list.dirs(file.path("../models_2015/sitka"),recursive=TRUE)

.getModelOutput <- function(fp) {
	#filter on the correlation file to ensure model converged
	if(file.exists(paste0(fp,"/ham.cor"))) {
		output <- read.admb(paste0(fp,"/ham"))
		output$Model <- fp
		return (output)
	} else {
		return (NULL)
	}
}

M <- lapply(.MODELDIRS,.getModelOutput)
M <- M[!sapply(M, is.null)]  # Remove NULL models from list (unconverged).


# ---------------------------------------------------------------------------- #
# PLOT DATA INPUTS
# ---------------------------------------------------------------------------- #

plot.all.data <- function(D) {
	# print(names(D))
	p.data <- NULL
	p.data$catch   <- plot.catch(D, nm = "data_ct_raw", y="Catch (tons)")
	p.data$egg <- plot.eggIndex(D,fit=FALSE) 
	p.data$log_egg <- plot.eggIndex(D,fit=FALSE,log.scale=TRUE) 
	p.data$sp_waa  <- plot.waa(D, nm = "data_sp_waa")
	p.data$cm_waa  <- plot.waa(D, nm = "data_cm_waa")
	p.data$sp_comp <- plot.comp(D, nm = "data_sp_comp")
	p.data$cm_comp <- plot.comp(D, nm = "data_cm_comp")

	cat("Finished plot.all.data\n")
	return( p.data )
}

# P.data <- plot.all.data(M[[1]])

data.plots <- lapply(M,plot.all.data)


# ---------------------------------------------------------------------------- #
# PLOT MLE MODEL OUTPUT
# ---------------------------------------------------------------------------- #

plot.all.output <- function(D) {
	p.output <- NULL
	p.output$label <- D$Model
	p.output$ssb <- plot.ssb(D)
	p.output$ft  <- plot.ft(D)
	p.output$ctfit<- plot.eggIndex(D,sfx="catch",fit=TRUE)
	p.output$egg <- plot.eggIndex(D,fit=TRUE)  
	p.output$logegg <- plot.eggIndex(D,fit=TRUE,log.scale=TRUE)  

	p.output$resd_cct <- plot.resd(D,nm="resd_catch")
	p.output$resd_egg <- plot.resd(D,nm="resd_egg_dep")
	p.output$resd_rec <- plot.resd(D,nm="resd_rec")
	p.output$resd_csp <- plot.resd(D,nm="resd_sp_comp")
	p.output$resd_ccm <- plot.resd(D,nm="resd_cm_comp")

	cat("Finished plot.all.output\n")
	return(p.output)
}

output.plots <- lapply(M,plot.all.output)

mle <- output.plots[[1]]
lapply(names(mle[-1]), 
  function(x)ggsave(filename=paste("../docs/modelDescription/figs/",x,".png",sep=""), plot=mle[[x]]))

# ---------------------------------------------------------------------------- #
# PLOT SIMULATION RESULTS
# ---------------------------------------------------------------------------- #

plot.all.simulations <- function(D) {

	# get index for simulation directories
	id <- grep("sims",sapply(D,"[[","Model"))

	# select only simulation model runs
	T  <- D[-id]
	S  <- D[id]


	p.sims <- NULL
	p.sims$label <- D$Model


	# Spawning stock biomass 
	yrs <- D[[1]]$year
	ssb <- cbind(Year=yrs,as.data.frame(sapply(S,"[[","ssb"))) %>% 
				 gather(Simulation,SSB,-Year)
	sitka <- cbind(Year=yrs,as.data.frame(sapply(T,"[[","ssb"))) %>%
				 gather(Simulation,SSB,-Year)

	print(head(ssb))
	# p.sims$ssb <- ggplot(ssb,aes(Year,SSB,color=factor(Simulation))) + 
	# 							geom_line(alpha=0.8) + 
	# 							geom_line(aes(Year,SSB),alpha=0.3,size=1.5,data=sitka,color="black") + 
	# 							labs(x="Year",y="Spawning stock biomass (mt)",color="Simulation")

	p.sims$ssb <- plot.ssb(T[[1]]) + 
								geom_line(aes(Year,SSB,color=factor(Simulation)),alpha=0.5,size=0.5,data=ssb) + 
								labs(x="Year",y="Spawning stock biomass (mt)",color="Simulation")


	# Fishing mortality rates
	fft <- cbind(Year=yrs,as.data.frame(sapply(S,"[[","ft"))) %>%
				 gather(Simulation,fft,-Year)
	# S
	p.sims$ft <- plot.ft(T[[1]]) +
							 geom_line(data=fft,aes(Year,fft,color=factor(Simulation)),alpha=0.5,size=0.5)+ 
							 labs(x="Year",y="Instantaneous fishing mortality rate")

  # Box plots for theta
	nm <- c("M","Rinit","Rbar","Ro","Reck","SigmaR")
	theta.pred <- sapply(D,"[[","theta")
	theta.true <- sapply(D,"[[","theta_ival")
	theta.ln2r <- t(log(theta.true/theta.pred,base=2)[,-1])
	colnames(theta.ln2r) <- nm

	p.sims$theta <- boxplot(theta.ln2r,ylim=c(-1,1))

	return(p.sims)
}
	
simulation.plots <- plot.all.simulations(M)

