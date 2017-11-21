#Plot observed and predicted data.

plot.eggIndex <- function(D, sfx="egg_dep",fit=FALSE,log.scale=FALSE) { 
	
	if(sfx=="egg_dep") {
		ylbl <-  "Egg Deposition log(trillions)"
	} else if (sfx=="catch") {
		ylbl <-  "Catch (mt)"
	}

	# Observed data
	data <- paste0("data_",sfx)
	data <- as.data.frame(D[[data]])
	colnames(data) <- c("year","index","log.se")

	# Predicted observations
	if(fit) {
		pred <- paste0("pred_",sfx)
		pred <- as.data.frame(cbind(D$year,D[[pred]]))
		colnames(pred) <- c("year","pred")
	}


	# Add LCI and UCI to data.frame
	z  <- 1.96

	df <- data %>%
				transform(index = ifelse(index<=0,NA,index)) %>%
				mutate(ln.It = log(index)) %>%
				mutate(lci = exp(ln.It - z*log.se),
				       uci = exp(ln.It + z*log.se)) %>%
				tbl_df()

	if(fit) {
		df <- df %>% left_join(pred,by="year") %>% tbl_df()
	}


	p <- ggplot(df,aes(year,index)) +
	geom_pointrange(aes(ymin=lci,ymax=uci),size=0.5,fatten=1) +
	labs(x="Year",y=ylbl) +
	ggtitle(D$Model)

	if(fit) {
		p <- p + geom_line(aes(year,pred),color="blue")
	}

	if(log.scale) {
		p <- p + scale_y_log10() + labs(y= ylbl )
	}
	return(p)
}