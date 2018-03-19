### Find 'string' in 'file' and read vector from next line
readVec <- function(string, file){
  txt <- readLines(file)
  skip <- match(string, txt)
  vec <- scan(file, quiet=TRUE, skip=skip, nlines=1)
  return(vec)
}

### Find 'string' in 'file' and read matrix with 'nrow' rows from next line
readMat <- function(string, file, nrow){
  txt <- readLines(file)
  skip <- match(string, txt)
  mat <- as.matrix(read.table(file, skip=skip, nrows=nrow))
  dimnames(mat) <- NULL
  return(mat)
}

## Function to read a basic AD Model Builder fit.
## Use for instance by:
## simple.fit <- readFit('c:/admb/examples/simple')
## Then the object 'simple.fit' is a list containing sub-objects
# 'names', 'est', 'std', 'cor', and 'cov' for all model
# parameters and sdreport quantities.
readFit <- function(file){
	ret <- list()
	parfile <- as.numeric(scan(paste(file,'.par', sep=''), what='', n=16, quiet=TRUE)[c(6,11,16)])
	ret$nopar <- as.integer(parfile[1])
	ret$nlogl <- parfile[2]
	ret$maxgrad <- parfile[3]
	file <- paste(file,'.cor', sep='')
	lin <- readLines(file)
	ret$npar <- length(lin)-2
	ret$logDetHess <- as.numeric(strsplit(lin[1], '=')[[1]][2])
	sublin <- lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!=''])
	ret$names <- unlist(lapply(sublin, function(x) x[2]))
	ret$est <- as.numeric(unlist(lapply(sublin, function(x) x[3])))
	ret$std <- as.numeric(unlist(lapply(sublin, function(x) x[4])))
	ret$cor <- matrix(NA, ret$npar, ret$npar)
	corvec <- unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)]))
	ret$cor[upper.tri(ret$cor, diag=TRUE)] <- as.numeric(corvec)
	ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)]
	ret$cov <- ret$cor*(ret$std%o%ret$std)
	return(ret)
}


## function used to calculate equilibrium depletion per one recruit
## can be used with uniroot to calculate the F that would result in equilibrium depletion
calc_equil_depl <- function(ages, W_a, M, F, ref=FALSE){

	Naf <- calc_equil_numbers(ages=ages, M=M, F=F, R0=1)
	Na0 <- calc_equil_numbers(ages=ages, M=M, F=0, R0=1)

        Nofish <- sum(Na0*W_a)
        Fish <- sum(Naf*W_a)

        ## automatically returns SPR
        ratio <- Fish/Nofish
        if(ref==FALSE) return(ratio)
            
        ## can use uniroot function on call to calc_ref to calculate the fishing mortality rate that results in a specified ratio, then compare current fishing mortality with this reference point
        if(ref!=FALSE){
            diff <- ref - ratio
            return(diff)
        }
}

calc_equil_numbers <- function(ages, M, F, R0){
	Na <- rep(NA, length(ages))
	Na[1] <- R0

	for(i in 2:length(ages)){
		if(i<length(ages)) Na[i] <- Na[i-1]*exp(-M-F)
		if(i==length(ages)) Na[i] <- (Na[i-1]*exp(-M-F))/(1-exp(-M-F))
	}
	return(Na)
}

simPop <- function(linf, 
					vbk, 
					t0, 
					lwa, 
					lwb, 
					M, 
					AgeMax=NULL, 
					S50,
					S95,
					M50,
					M95,
					R0,
					h,
					SigmaR,
					SigmaF,
					Nyear,
					InitialDepl,
					PropFcrash,
					EffortDyn,
					seed){

	##################
	## Age and growth
	##################
	if(is.null(AgeMax)) AgeMax <- ceiling(-log(0.01)/M)
	ages <- 0:AgeMax

	L_a <- linf * (1 - exp(-vbk*(ages - t0)))
	W_a <- lwa * L_a ^ lwb

	###################
	## Process error
	###################
	set.seed(seed)
	RecDev <- rnorm(Nyear, mean=-(SigmaR^2)/2, sd=SigmaR)
	EffortDev <- rnorm(Nyear, mean=-(SigmaF ^ 2)/2, sd=SigmaF)

	####################
	## Reference points
	####################
	Finit <- uniroot(calc_equil_depl,
						lower=0,
						upper=5,
						ages=ages,
						W_a=W_a,
						M=M,
						ref=InitialDepl)$root

	Fcrash <- uniroot(calc_equil_depl,
						lower=0,
						upper=10,
						ages=ages,
						W_a=W_a,
						M=M,
						ref=0.05)$root

	Fterm <- Fcrash * PropFcrash

	####################
	## Effort dynamics
	####################

	if(EffortDyn == "Constant"){
		E_t <- rep(1, Nyear)
		q <- Finit/max(E_t)
	}
	if(EffortDyn == "One-way"){
		E_t <- seq(0.01,by=0.05,length=Nyear)
		q <- Fterm/max(E_t)
	}
	if(EffortDyn == "Two-way"){
		E1_yr <- ceiling(Nyear/3)
		E2_yr <- ceiling(Nyear/3)
		E3_yr <- Nyear - E1_yr - E2_yr
		E_t1 <- seq(0.01, by=0.05,length=E1_yr)
		E_t2 <- rep(E_t1[length(E_t1)], E2_yr)
		E_t3 <- seq(from=E_t1[length(E_t1)], to=E_t1[length(E_t1)]/2, length=E3_yr)
		E_t <- c(E_t1, E_t2, E_t3)
		q <- Fterm/max(E_t)
	}
	F_t <- q * E_t * exp(EffortDev)

	###########################
	## Selectivity & Maturity
	###########################
	S_a <- 1 / (1 + exp(-log(19) * (ages - S50) / (S95 - S50)))
	Mat_a <- 1 / (1 + exp(-log(19) * (ages - M50) / (M95 - M50)))

	####################
	## Recruitment
	####################
	R_t <- rep(NA, Nyear)
	R0_inp <- R0 * exp(RecDev[1])
	R_t[1] <- R0_inp

	####################
	## Initialization
	####################
	N_at <- F_at <- matrix(NA, nrow=length(ages), ncol=Nyear)
	SB_t <- VB_t <- rep(NA, Nyear)

	####################
	## Mortality
	####################
	for(y in 1:Nyear){
		for(a in 1:length(ages)){
			F_at[a,y] <- F_t[y] * S_a[a]
		}
	}
	Z_at <- F_at + M

	####################
	## Population
	####################

	## initialize equilibrium numbers
	N_at[,1] <- calc_equil_numbers(ages=ages, M=M, F=F_t[1], R0=R_t[1])

	## initialize spawning biomass
	SB0 <- sum(R0 * exp(-M * ages) * Mat_a * W_a)
	SB_t[1] <- sum(N_at[,1] * Mat_a * W_a)
	VB_t[1] <- sum(N_at[,1] * S_a * W_a)

	## project the numbers matrix
	for(y in 2:Nyear){
		R_t[y] <- ((4 * h * R0 * SB_t[y-1]) / (SB0 * (1-h) + SB_t[y-1] * (5*h - 1))) * exp(RecDev[y])
		for(a in 1:length(ages)){
			if(a==1) N_at[a,y] <- R_t[y]
			if(a>1 & a<length(ages)) N_at[a,y] <- N_at[a-1,y-1] * exp(-Z_at[a-1,y-1])
			if(a==length(ages)) N_at[a,y] <- N_at[a-1,y-1] * exp(-Z_at[a-1,y-1]) + N_at[a,y-1] * exp(-Z_at[a,y-1])	
		}
		SB_t[y] <- sum(N_at[,y] * Mat_a * W_a)
		VB_t[y] <- sum(N_at[,y] * S_a * W_a)
	}
	N_t <- colSums(N_at[-1,])

	####################
	## Observation model
	####################

	C_at <- matrix(NA, nrow=length(ages), ncol=Nyear)
	## catch in weight
	for(y in 1:Nyear){
		for(a in 1:length(ages)){
			C_at[a,y] <- N_at[a,y] * W_a[a] * (F_at[a,y] / Z_at[a,y]) * (1 - exp(-Z_at[a,y])) 
		}
	}
	C_t <- colSums(C_at)

	CPUE_t <- C_t / E_t

	page <- matrix(NA, nrow=length(ages), ncol=Nyear)
	for(y in 1:Nyear){
		page[,y] <- N_at[,y] * S_a
	}
	page <- sapply(1:Nyear, function(x) page[,x]/sum(page[,x]))

	####################
	## Return output
	####################

	out <- NULL
	out$C_t <- C_t
	out$CPUE_t <- CPUE_t
	out$page <- page
	out$N_t <- N_t
	out$VB_t <- VB_t
	out$SB_t <- SB_t
	out$SB0 <- SB0
	out$R_t <- R_t
	out$E_t <- E_t
	out$F_t <- F_t
	out$q <- q
	out$S_a <- S_a
	out$Mat_a <- Mat_a
	out$L_a <- L_a
	out$W_a <- W_a
	out$ages <- ages
	out$Finit <- Finit
	out$Fterm <- Fterm
	out$Fcrash <- Fcrash
	out$RecDev <- RecDev
	out$EffortDev <- EffortDev
	out$Nyear <- Nyear
	out$Nage <- length(ages)
	out$M <- M
	out$S50 <- S50
	out$S95 <- S95
	out$M50 <- M50
	out$M95 <- M95
	out$SigmaR <- SigmaR
	out$SigmaF <- SigmaF
	out$R0 <- R0
	out$h <- h
	return(out)
}