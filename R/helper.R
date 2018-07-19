
# Libraries, ggplot themes, and user-defined fxns.
# Author: Jane Sullivan
# Contact: jane.sullivan1@alaska.gov

# libraries ----

if(!require("mosaic"))   install.packages("mosaic") # derivedFactor, derivedVariable. masks over a lot of fxns, but generally only improves their utility
if(!require("tidyverse"))   install.packages("tidyverse") # dplyr, ggplot, etc.
if(!require("lubridate"))   install.packages("lubridate") # dates functions like yday, dmy/mdy
# if(!require("mgcv"))   install.packages("mgcv") # gams
# if(!require("gridExtra"))   install.packages("gridExtra") # multipanneled plots
if(!require("data.table"))   install.packages("data.table") # dcast, foverlaps
# if(!require("ROracle"))   install.packages("ROracle") # database access through R
# if(!require("broom"))   install.packages("broom") # tidying regression model output
# if(!require("padr"))   install.packages("padr") # fills in missing values in a time series
if(!require("knitr"))   install.packages("knitr") # r markdown
# if(!require("forcats"))   install.packages("forcats") # releveling factors
if(!require("cowplot"))   install.packages("cowplot") # plot_grid and so much else
if(!require("R2admb"))   install.packages("R2admb") # run admb from r
if(!require("ggthemes"))   install.packages("ggthemes") # access to 150 colour palettes from canva.com design school, among other things
# https://www.canva.com/learn/100-color-combinations/
# http://makeadifferencewithdata.com/wp-content/uploads/2016/12/color-palettes.txt
# if(!require("scales"))   install.packages("scales") # used to expand colour palettes
# install.packages("devtools")
# devtools::install_github("ben-williams/FNGr")
# library("FNGr")

# ggplot themes ----

windowsFonts(Times=windowsFont("Times New Roman"))

theme_sleek <- function(base_size = 12, base_family = "Times") {
  half_line <- base_size/2
  theme_light(base_size = 12, base_family = "Times") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      #axis.text = element_text(colour = "grey30"),
      #axis.title = element_text(colour = "grey30"),
      #legend.title = element_text(colour = "grey30"),#, size = rel(0.9)
      panel.border = element_rect(fill = NA),#, colour = "grey70", size = 1),
      legend.key.size = unit(0.9, "lines"),
      #legend.text = element_text(size = rel(0.7)),#, colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)#,
      #plot.title = element_text(colour = "grey30"),#, size = rel(1)
      #plot.subtitle = element_text(colour = "grey30")#, size = rel(.85)
    )
}

theme_set(theme_sleek())

# User-defined functions ----

# tickr() - makes nice ticks and breaks on your ggplots

# Depends on dplyr
tickr <- function(
  data, # dataframe
  var, # column of interest
  to # break point definition 
){
  
  VAR <- enquo(var) # makes VAR a dynamic variable
  
  data %>% 
    distinct(!!VAR) %>%
    ungroup(!!VAR) %>% 
    mutate(labels = ifelse(!!VAR %in% seq(to * round(min(!!VAR) / to), max(!!VAR), to),
                           !!VAR, "")) %>%
    select(breaks = UQ(VAR), labels)
}

# Functions from SJD Martell, formerly in globals.r

read_admb <-
  function(ifile)
  {	
    ret=read.fit(ifile)
    
    fn=paste(ifile,'.rep', sep='')
    A=read.rep(fn)
    A$fit=ret
    
    pfn=paste(ifile,'.psv',sep='')
    if(file.exists(pfn))
      A$post.samp=read.psv(pfn)
    
    return(A)
  }

read_fit <-
  function(ifile)
  {
    # __Example:             
    #	file <-("~/admb/simple")
    #	A <- reptoRlist(file)
    #	Note there is no extension on the file name.
    
    ## The following is a contribution from:
    ## Anders Nielsen that reads the par & cor files.
    ret<-list() 
    parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
                             what='', n=16, quiet=TRUE)[c(6,11,16)]) 
    ret$nopar<-as.integer(parfile[1]) 
    ret$nlogl<-parfile[2] 
    ret$maxgrad<-parfile[3] 
    file<-paste(ifile,'.cor', sep='') 
    lin<-readLines(file) 
    ret$npar<-length(lin)-2 
    ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
    sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
    ret$names<-unlist(lapply(sublin,function(x)x[2])) 
    ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
    ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
    ret$cor<-matrix(NA, ret$npar, ret$npar) 
    corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
    ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
    ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
    ret$cov<-ret$cor*(ret$std%o%ret$std)
    return(ret)
  }

read_rep <- 
  function(fn)
  {
    # The following reads a report file
    # Then the 'A' object contains a list structure
    # with all the elemements in the report file.
    # In the REPORT_SECTION of the AMDB template use 
    # the following format to output objects:
    #  	report<<"object \n"<<object<<endl;
    #
    # The part in quotations becomes the list name.
    # Created By Steven Martell
    options(warn=-1)  #Suppress the NA message in the coercion to double
    
    
    ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
    idx=sapply(as.double(ifile),is.na)
    vnam=ifile[idx] #list names
    nv=length(vnam) #number of objects
    A=list()
    ir=0
    for(i in 1:nv)
    {
      ir=match(vnam[i],ifile)
      if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
      dum=NA
      if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
      if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))
      
      if(is.numeric(dum))#Logical test to ensure dealing with numbers
      {
        A[[vnam[i]]]=dum
      }
    }
    options(warn=0)
    
    return(A)
  }

read_psv <-
  function(fn, nsamples=10000)
  {
    #This function reads the binary output from ADMB
    #-mcsave command line option.
    #fn = paste(ifile,'.psv',sep='')
    filen <- file(fn, "rb")
    nopar <- readBin(filen, what = integer(), n = 1)
    mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
    mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
    close(filen)
    return(mcmc)
  }

# A simple function for creating transparent colors
# Author: Nathan Stephens (hacks package)
colr <- 
  function(col.pal=1,a=1)
  {
    col.rgb<-col2rgb(col.pal)/255
    rgb(t(col.rgb),alpha=a)
  }

