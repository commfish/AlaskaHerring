# Rscript for your model outputs.

.REP <- "ham"
# setwd("/Users/stevenmartell1/Library/Mobile Documents/com~apple~CloudDocs/icSEASTATE/ADFG_HerringModel/models_2015/sitka")
source("../../R/Globals.R")
# source("../../R/Ham.R")

if(!file.exists("ham.par")) {
	# system("ham.exe -ind ham.dat")  #windoz box
	system("make")
}
S <- read.admb("ham")


# plot.datafit(S, sfx="egg_dep", fit=TRUE)