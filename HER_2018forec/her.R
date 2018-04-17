
# HER - ADFG's herring model. Original ADMB code written by SJD Martell. Helper
# files and documention contributed by M. Rudd.

source("R/tools.R") 
library(R2admb)


# clean up directory - remove unnecessary files to run
# add ".pin" when wanting to save .pin file, or any other file to save
setwd("HER_2018forec/")
need_files <- c(".tpl", ".ctl", ".dat", ".R")
files_present <- list.files()
keep_index <- unlist(sapply(1:length(need_files), function(x) which(grepl(need_files[x], files_present))))
rm_files <- files_present[-keep_index]
remove <- sapply(1:length(rm_files), function(x) unlink(rm_files[x], TRUE))
# setwd("..") # backs out one folder

## compile model
setup_admb()
compile_admb("her", verbose=TRUE)

## run MAP
run_admb("her")

years <- 1971:2015

#Maturity
readMat("mat", file="her.rep", nrow = length(years))
readVec("mat_params[1]", file="her.rep")
readVec("mat_params[2]", file="her.rep")

# Natural mortality
natmat <- readMat("Mij", file = "her.rep", nrow = length(years))[,1]
plot(natmat ~ years, ylim=c(0, max(natmat)*1.1), type="l", lwd=2)

## read report from initial MAP run - maximum a posteriori estimation (i.e.
## maximum likelihood using priors!)
ssb <- readVec("ssb", file="her.rep")
## put in thousands
ssb <- ssb/1000

plot(ssb ~ years, ylim=c(0, max(ssb)*1.1), type="l", lwd=2)
## run simulation with seed 123
run_admb("her", extra.args="-sim 123")

## run MCMC
run_admb("her", extra.args="-mcmc 10000 -mcsave 10")
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
