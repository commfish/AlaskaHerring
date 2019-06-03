# simulation.R
library(dplyr)
library(tidyr)
library(ggplot2)
# write the code in here to read all 00001. directories and sumarize results.
# setwd("/Users/steve/Library/Mobile Documents/com~apple~CloudDocs/icSEASTATE/ADFG_HerringModel/app")
source(file.path("./globals.R"))

readOutput <- function(d) {
	if(file.exists(file.path(d,"ham.cor"))){
		return(read.admb(file.path(d,"ham")));
	}
	else {
		return(NULL)
	}
}

dn <- dir(path="../models_2015/sitka/sims",pattern="^[[:digit:]]",full.names=TRUE);
nf <- paste0("sims",".Rdata");
sims <- lapply(dn,readOutput)
save(sims,file=nf)


# Create Data Frame of Parameter Estimates.
t1 <- lapply(sims,'[[','fit')
t2 <- lapply(t1,'[[','est')
df <- do.call(rbind,t2)

theta <- sims[[1]]$theta_ival

# boxplot(t(log2(t(df[,1:5])/theta)),ylim=c(-1,1))

# A nicer boxplot.
str   <- c("M","Ri","Rbar","Ro","Kappa")
ndf   <- df[,1:5]; colnames(ndf) <- str
err   <- as.data.frame(t((t(ndf)-theta)/theta)*100)

p <- ggplot(gather(err),aes(key,value)) + geom_boxplot() + 
			ylim(c(-100,100)) +
		labs(x="Parameter",y="Relative Error (%)")