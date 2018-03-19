main <- "C:\\merrill\\ADMB_course" 
source(file.path(main, "tools.R"))

wd <- "C:\\merrill\\ADMB_course\\HER"
setwd(wd)

## clean up directory - remove unnecessary files to run
## add ".pin" when wanting to save .pin file, or any other file to save
need_files <- c(".tpl", ".ctl", ".dat", ".R") 
files_present <- list.files(wd)
keep_index <- unlist(sapply(1:length(need_files), function(x) which(grepl(need_files[x], files_present))))
rm_files <- files_present[-keep_index]
remove <- sapply(1:length(rm_files), function(x) unlink(rm_files[x], TRUE))

library(R2admb)

## compile model
compile_admb("her", verbose=TRUE)

## run MAP
run_admb("her")

## read report from initial MAP run
ssb <- readVec("ssb", file="her.rep")
plot(ssb, ylim=c(0, max(ssb)*1.1), type="l", lwd=2)

## run simulation with seed 123
run_admb("her", extra.args="-sim 123")

## run MCMC
run_admb("her", extra.args="-mcmc 10000 -mcsave 10")
run_admb("her", extra.args="-mceval")

## posterior distributions
ssb_ps <- read.table(file.path(wd, "ssb.ps"))
natural_ps <- read.table(file.path(wd, "natural.ps"), header=TRUE)

par(mfrow=c(2,1))
plot(natural_ps[,1])
abline(h=median(natural_ps[,1]), col="red")
plot(natural_ps[,2])
abline(h=median(natural_ps[,2]), col="red")


## from tech doc:
## 1) first fit to sitka data
## 2) save .par file as her.pin
