# ---------------------------------------------------------------------------- #
# Herring Age Model, or ham 
# June 29, 2016
# Author: Steven Martell
# email: steve@seastateinc.com
# phone: 808 280-2548
# ---------------------------------------------------------------------------- #
library(ggplot2)
library(dplyr)
library(tidyr)

# Read in the data from the model report, par, and cor files.
source(file.path("./globals.R"))
# D <- read.admb("../models_2015/sitka/ham")
# sb.file <- "../models_2015/sitka/ssb.ps"
D <- read.admb("../models_2015/craig/ham")
sb.file <- "../models_2015/craig/ssb.ps"
if(file.exists(sb.file)){
	D$post.samp.ssb=read.table(sb.file)
	colnames(D$post.samp.ssb) <- paste0("year",D$year)
}

# ---------------------------------------------------------------------------- #
# DATA SECTION
# ---------------------------------------------------------------------------- #

plot.catch <- function(D=D, nm = "data_ct_raw",...) {
	df <- as.data.frame(D[[nm]])
	colnames(df) <- c("Year","Catch","log.se")
	z  <- 1.96
	df <- df %>% 
				mutate(ln.ct = log(Catch)) %>%
				mutate(lci = exp(ln.ct - z*log.se),
				       uci = exp(ln.ct + z*log.se)) %>%
				mutate(std = 1.96*sqrt(log(log.se+1))) %>%
				mutate(lower=Catch-std*Catch,upper=Catch+std*Catch)

	ggplot(df,aes(Year,Catch)) +
		geom_pointrange(aes(ymin = lci, ymax = uci),size=0.5,fatten=2) + 
		labs(x="Year",...) + ggtitle(D$Model)
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

# Deprecate
plot.datafit <- function(D=D, sfx="egg_dep", fit=FALSE) {
	data <- paste0("data_",sfx)
	data <- as.data.frame(D[[data]])
	colnames(data) <- c("year","index","log.se")

	pred <- paste0("pred_",sfx)
	pred <- as.data.frame(cbind(D$year,D[[pred]]))
	colnames(pred) <- c("year","pred")

	df   <- right_join(data,pred) %>%
					mutate(ln.index=log(index)) %>%
					mutate(std = 1.96*sqrt(log(log.se+1))) %>%
					mutate(lower=exp(ln.index-std*ln.index),
					       upper=exp(ln.index+std*ln.index)) %>%
					tbl_df()

	# print(head(df))

	ggplot(df,aes(year,index)) + 
	geom_pointrange(aes(ymin = lower, ymax = upper),size=0.5,fatten=2) +
	labs(x="Year",y="Egg Deposition (trillions)") +
	ggtitle(D$Model) +
	if(fit) geom_line(aes(year,pred),alpha=0.8) 

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
plot.ssb.post <- function(D=D) {
	# if(file.exists(sb.file)){
		
		ssb.ps <- D$post.samp.ssb
		colnames(ssb.ps) <- paste(D$year)
		df <- gather(ssb.ps,Year,SSB)
		ggplot(df,aes(Year,SSB/1000)) + 
		geom_violin(alpha=0.25,fill="steel blue",size=0.25,
	            draw_quantiles = c(0.25, 0.5, 0.75)) +
		ylim(c(0,NA)) + labs(x="Year",y="Spawning Stock Biomass (1000 t)") +
		scale_x_discrete(breaks=seq(D$mod_syr,D$mod_nyr,by=5))

	# }
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
d1 <- plot.catch(D, nm = "data_ct_raw",y="Catch (tons)")
d2 <- plot.waa(D,nm = "data_sp_waa")
d3 <- plot.waa(D,nm = "data_cm_waa")
d4 <- plot.comp(D,nm = "data_cm_comp")
d5 <- plot.comp(D,nm = "data_sp_comp")




