#Created and maintained by Sara Miller Alaska Department of Fish & Game (sara.miller@alaska.gov) (August 2015)
#NOTES: For each stock, the MAXy, axis text label size, and width and height of figure may need to be adjusted. 
#To change to a New stock:
#(1) Line 125; setwd("H:/Herring/ADMB to R Figures/Herring_Models_August 2015/Working Folder/Sitka") #set path to where ADMB puts model results
#(2) Line 127; A = reptoRlist(fn="H:/Herring/ADMB to R Figures/Herring_Models_August 2015/Working Folder/Sitka/model.rep") #change based on output of best model
#(3) Line 162 and 163; Figure 1: change the U and L based on the best model and std report output (only for Sitka and Craig). If running other stocks,
    #Line 176;"expr=errbar(C$Year,C$tot_obs_egg,C$tot_est_egg+U,C$tot_est_egg-L,cap=.01,add=TRUE,pch=1,cex=0,col="light grey",errbar.col="blue",lwd=1)"
    #is not needed
#(4) Line 602 and 649; In Figure 8, make sure the right stock's data is being used. Also, make sure the "facet_wrap(~Year,nrow=7, ncol=6)" states the correct row and cols.
#(5) Line 712 and 758; In Figure 9,make sure the right stock's data is being used. Also, make sure the "facet_wrap(~Year,nrow=7, ncol=6)" states the correct row and cols.
     #All other stock data should be commented out with a "#".
#(6) There is an alternate figure 3 for Tenakee since the past years are consecutive (lines 1067-1105). Run for Tenakee only.
#(7) Don't run the section for Sitka only (lines 1114-1270) unless you are running the Sitka stock. 
#(8) All occurrences of the stock name in Section III also need to be changed to the correct stock being run.
#SECTION #I:READ IN ADMB REPORT OUTPUT & SET UP FILES
#SECTION #II:CREATE FIGURES 
        #*Note: Mature is the same as pre-fishery and spawning is the same as post-fishery
        #Figure 1:  Survey- and model-estimated egg deposition.
        #Figure 2:  Survey-estimated spawning biomass plus catch (tons), model-estimated mature biomass (tons), and model-estimated mature biomass forecast (tons).
        #Figure 3:  Comparison of past and current survey-estimated mature biomass (survey-estimated spawning biomass plus catch), model- estimated mature biomass, and model- estimated mature biomass forecasts (tons).
        #Figure 4:  Residuals from model fits to survey egg deposition and Ricker spawner-recruit function.
        #Figure 5:  Model estimates of age-3 recruit strength (numbers of age-3 mature and immature fish).
        #Figure 6:  Spawning population biomass (blue bars;top figure), spawning population abundance (blue bars;middle figure), population abundance (immature and spawning abundance) (blue bars;bottom figure), and commercial fishery harvest (yellow bars) over time. The combination of the blue and yellow bars (total height of each bar) is the mature biomass, mature population abundance, or total population abundance.
        #Figure 7a-Figure 7c:  Model estimates of gear selectivity at age (a), maturity at age (b), and survival (c) by year.
        #Figure 8:  Observed seine or gillnet (red line with square points) and model-estimated (bar) catch-age composition.
        #Figure 9:  Observed cast net (red line with square points) and model-estimated (bar) spawning-age composition.
        #Figure 10: Spawning biomass (tons) versus age-3 abundance (millions of mature and immature fish) (blue circles) with Ricker-estimated age-3 abundance (red triangles), from 1991-2014. 
        #Figure 11: Projected mature biomass at age (tons) for forecast year.
        #Figure 12: Projected percentage of mature numbers at age for forecast year. 
        #Figure 13: Forecasted weight at age
        #Figure 14: Stacked bar graph of catch (orange), spawning biomass (green), GHL (blue), and the spawning biomass forecast (pink). The harvest (or GHL) plus the spawning biomass equals the mature biomass. 
        #Figure 15: Spawning age composition residuals.
        #Figure 16: Estimated spawning age composition.
        #Figure 17: Estimated catch age composition.
#EXTRA FIGURES FOR SITKA ONLY
        #Figure 8a (Sitka), Figure 9a (Sitka)
#SECTION #III: OUTPUT FILES TO ONE EXCEL WORKBOOK & FORMAT
#######################################################################################################################################
##LIBRARY DIRECTORY FOR DOWNLOADS - paste into library routing if needed to avoid duplicate "user library"
## - REQUIRES READ/WRITE ADMIN on lib="C:\~\R\win-library\3.1
##lib="C:\~\R\win-library\3.1"
###LIBRARIES -these are standard libraries; not all of them are always used
rm(list=ls(all=T))#Remove previous variables.
#install.packages("plyr")
#install.packages("reshape2")
#install.packages("lattice")
#install.packages("latticeExtra")
#install.packages("gridExtra")
#install.packages("ggplot2")
#install.packages("Hmisc")
#install.packages("MASS")
#install.packages("survival")
#install.packages("scatterplot3d")
#install.packages("vcd")
#install.packages("graphics")
#install.packages("calibrate")
#install.packages("scales")
#install.packages("xlsx")
#install.packages("extrafont")
library(plyr)
library(reshape2)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(Hmisc)
library(MASS)
library(survival)
library(scatterplot3d)
library(vcd)
library(grid)
library(calibrate)
library(scales)
library(extrafont)
library(xlsx)
library(RColorBrewer)
#font_import() #only do this one time - it takes a while
loadfonts(device="pdf")#http://gettinggeneticsdone.blogspot.com/2010/03/arrange-multiple-ggplot2-plots-in-same.html (code to put 2 plot together)
require(grid) #the code below makes it possible to create a two-panel figure in ggplot
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)  
  n <- length(dots)	
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}	
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}	
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}        
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )  
  ii.p <- 1	
  for(ii.row in seq(1, nrow)){	
    ii.table.row <- ii.row		
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}		
    for(ii.col in seq(1, ncol)){			
      ii.table <- ii.p			
      if(ii.p > n) break			
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))			
      ii.p <- ii.p + 1}}}
#function to extract data from ADMB report
reptoRlist = function(fn)
{
  ifile=scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx]  #list names
  nv=length(vnam)           #number of objects
  A=list()
  ir=0
  for(i in 1:nv)
  {
    ir=match(vnam[i],ifile)
    if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
    dum=NA
    if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
    if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
    
    if(is.numeric(dum))#Logical test to ensure dealing with numbers
    {
      A[[ vnam[i ] ]]=dum
    }
  }
  return(A)
}
#formatting axis function
fmt <- function(){
  function(x) format(x,nsmall = 2,scientific = FALSE)}
round2<-function(x){trunc(x+0.5)} 
##SET WORKING DIRECTORY-NEEDS TO BE CHANGED!!!!
# setwd("C:/2016 GIT/Sitka/ADMB/211a") #set path to where ADMB puts model results
path<-getwd()
A = reptoRlist(fn="model.rep") #change based on output of best model
# source("C:/2016 GIT/Sitka/ADMB/211a/sit.gz")#change based on output of best model (this file is needed for Figures 16 & 17)
#####################################################################################################################################
#####################################################################################################################################
#SECTION #I:READ IN ADMB REPORT OUTPUT & SET UP FILES; 3 csv files are created
#####################################################################################################################################
#####################################################################################################################################
##READ IN DATA FROM ADMB
#(1) FIGDATA
FIGDATA<- read.csv("FIGDATA.dat", header=TRUE, sep="") 
write.csv(FIGDATA, "FIGDATA.csv") 
FIGDATA<- read.csv("FIGDATA.csv", header=TRUE, stringsAsFactors = FALSE) #Calculate certain  variables for figures
FIGDATA[FIGDATA==-9] <- NA
FIGDATA$tcb[is.na(FIGDATA$tcb)] <- 0 

FIGDATA["tot_mat_B_tons"]<-FIGDATA$tot_mat_B/0.90718 #convert tonnes to tons
FIGDATA["tot_sp_B_tons"]<-FIGDATA$tot_sp_B/0.90718 #convert tonnes to tons

#(2) FIGDATAAGE
FIGDATAAGE<- read.csv("FIGDATAAGE.dat", header=TRUE, sep="") 
write.csv(FIGDATAAGE, "FIGDATAAGE.csv") 
FIGDATAAGE<- read.csv("FIGDATAAGE.csv", header=TRUE) 
FIGDATAAGE["Age"] <- c(3,4,5,6,7,8)#ages 3-8
FIGDATAAGE["Age2"] <- ifelse(FIGDATAAGE$Age>=8,"8+",FIGDATAAGE$Age) #add ages 8+
FIGDATAAGE["for_mat_baa_tons"]<-FIGDATAAGE$for_mat_baa/0.90718 #convert to tons
###################################################################################################################################
###################################################################################################################################
#SECTION #II:CREATE FIGURES 
###################################################################################################################################

###################################################################################################################################
#GRAPHIC ONE: 
#Figure 1: Survey- and model-estimated egg deposition.
C<-subset(FIGDATA, select=c(Year, tot_obs_egg,tot_est_egg)) 
library(Hmisc)
#error bars for Stika and Craig only (from the egg variance sheet-upper and lower CI)
U<-c(1.23,  1.17,	1.15,	1.00,	1.28,	1.18,	0.97,	1.74,	2.59,	1.03,	0.84,	1.259686251,	1.851147636,	1.751126689,	0.560987576,
1.508697833,	1.633749193,	1.695401525,	1.255367509,	2.387620843,	2.795968238,	2.215761696,	1.462234716,	2.289501604,
2.650921062,	6.384923885,	2.279693146,	2.872760889,	29.05686308,	3.863145759,	4.816134565,	4.222205391,	1.634805164,
4.043867944,	1.288746439,  1.332721825) #Fix 2015 values
L<-c(1.23,  1.17,	1.15,	1.00,	1.28,	1.18,	0.97,	1.74,	2.59,	1.03,	0.84,	0.985282715,	1.639191663,	1.382705136,	0.48874792,
1.398489625,1.265029243,	1.286735024,	1.146877561,	1.827147032,	2.534746454,	1.882753246,	1.475607608,	1.863883108,
2.277982827,	3.540565615,	1.707320794,	2.568958439,	14.54490887,	3.237073047,	3.7630828,	3.942674884,	1.578639917,
2.996229014,	1.089460882, 1.154768444) #Fix 2015 values

windowsFonts(A = windowsFont("Times New Roman"))
MAXy<-max(C$tot_obs_egg, C$tot_est_egg,na.rm=TRUE)*1.1 # set y max. limit
png(file='Figure 1.png', res=200, width=8, height=4, units ="in")  
op <- par(family = "Times New Roman")
plot(C$Year,C$tot_obs_egg,pch=21,col="black",xaxt="n", bg="blue",cex=1.2,lwd=1,ylab="Eggs spawned (trillions)",xlab="Year",
     cex.axis=1.2,cex.lab=1.2, ylim=c(0, MAXy), xaxt="n", family="A")
lines(C$Year,C$tot_obs_egg,lty=2,lwd=1,col="black") #observed
lines(C$Year,C$tot_est_egg,lwd=3,col="black") #predicted
points(C$Year,C$tot_est_egg,pch=17,col="red", cex=1.3) #predicted
points(C$Year,C$tot_obs_egg,pch=21,col="black", bg="blue", cex=1.3) #observed
axis(side=1,at=seq(min(C$Year),max(C$Year),1),cex.axis=1, las=2)

expr=errbar(C$Year,C$tot_obs_egg, C$tot_obs_egg+U,C$tot_obs_egg-L,cap=.01,add=TRUE,pch=1,cex=0,col="light grey",errbar.col="blue",lwd=1)

legend("topleft",c("Survey-estimated egg deposition","Model-estimated egg deposition"),pch=c(16,17),col=c("blue", "red"),
      cex=1, bty="n")

dev.off()
Figure_1_data<-C
rm (C, MAXy, U, L)
ls(all.names = TRUE)
#****************************************************************************************************************************************************
#GRAPHIC TWO: 
#Figure 2: Survey-estimated spawning biomass plus catch (tons), model-estimated mature biomass (tons), 
#and model-estimated mature biomass forecast (tons).
#ONLY THE MODEL YEARS SHOULD BE USED
num<-c(35000,30000,29500,23500,38500,31000,25000,46000,58500,27000,23000,23500,43351,37150,14941,34990,40827,28611,34942,
       44554,57988,58756,40366,55769,69907,101305,66111,84501,247088,110946,126230,161904,62518,103267,48561, 58183)

#update survey-estimated spawning biomass for each stock based on the current spawn deposition file '20[XX] Corrected escapement (tons) using 10% egg loss' column AF

FIGDATA["spawn_dep"]<-num #survey-estimated spawning biomass in tons

y<-A$for_mat_B_st#forecast mature biomass (tons)
x<-subset(FIGDATA, select=c(Year, tot_mat_B_tons))#estimated mature biomass
Y<-c(max(x$Year)+1)
B<- data.frame(Y)  
B["Year"]<-B$Y
B["tot_mat_B_tons"]<-0
B<-subset(B, select=c(Year, tot_mat_B_tons))
C<-rbind(x,B) #to create x axis with a greater x value than max(Year)

png(file='Figure 2.png', res=200, width=10, height=5, units ="in")
MAXy<-max(C$tot_mat_B_tons,FIGDATA$spawn_dep+FIGDATA$tcb,na.rm=TRUE)*1.2# set y max. limit
windowsFonts(A = windowsFont("Times New Roman"))
op <- par(family = "Times New Roman")
plot(C$Year,C$tot_mat_B_tons, pch=21,col="white", bg="white",ylab="Biomass (tons)", xlab="Year",
     cex.axis=1.2,cex.lab=1.2,xaxt="n", ylim=c(0,MAXy), family="A") #placeholder
points(FIGDATA$Year,FIGDATA$spawn_dep+FIGDATA$tcb,pch=21, col="black",bg="blue", cex=1.3) #based on model years
lines(FIGDATA$Year,FIGDATA$spawn_dep+FIGDATA$tcb,lty=2,lwd=1, col="black")
lines(FIGDATA$Year,FIGDATA$tot_mat_B_tons,lwd=3)
points(FIGDATA$Year,FIGDATA$spawn_dep+FIGDATA$tcb,pch=21, col="black",bg="blue", cex=1.3)#survey estimated mature biomass
points(FIGDATA$Year,FIGDATA$tot_mat_B_tons,pch=17, col="red", cex=1.3)
points(max(FIGDATA$Year)+1,y,pch=8,col="red",cex=1)
textxy(max(FIGDATA$Year)+1, y,round2(y), pos=3, cex=1)
lines(FIGDATA$Year,FIGDATA$Threshold,col="grey",lty=1, lwd=2)
axis(side=1,at=seq(min(C$Year),max(C$Year)+1,1),cex.axis=1, las=2)
legend("topleft",c("Survey-estimated spawning biomass+catch","Model-estimated mature biomass",
                    "Mature biomass forecast","Threshold"), lty=c(NA,NA,NA,1),pch=c(16,17,8,NA),
       lwd=c(NA,NA,NA,2), col=c("blue","red","red", "grey"), cex=1, bty="n")
dev.off()
C<-subset(FIGDATA, select=c(Year, tot_mat_B_tons, spawn_dep, tcb))
Figure_2_data<-C
rm (B,C,y,x,Y, MAXy)
ls(all.names = TRUE)
#***********************************************************************************************************************
#GRAPHIC THREE: 
#Figure 3: Comparison of past and current survey-estimated mature biomass (survey-estimated spawning biomass plus catch), model- estimated mature biomass, 
#and model- estimated mature biomass forecasts (tons). Tenakee does not have a consistent three year forecast.
num<-c(35000,30000,29500,23500,38500,31000,25000,46000,58500,27000,23000,23500,43351,37150,14941,34990,40827,28611,34942,
       44554,57988,58756,40366,55769,69907,101305,66111,84501,247088,110946,126230,161904,62518,103267,48561, 58183)#update for each stock

FIGDATA["spawn_dep"]<-num #survey-estimated spawning biomass in tons
png(file='Figure 3.png', res=200, width=11, height=5, units ="in")  
MAXy<-max(A$yminusoneFOR, A$yminusthreeFOR, A$yminustwoFOR, A$for_mat_B_st,FIGDATA$spawn_dep+FIGDATA$tcb,na.rm=TRUE)*1.2# set y max. limit
par(xpd=T, mar=par()$mar+c(0,0,0,0))
windowsFonts(A = windowsFont("Times New Roman"))
op <- par(family = "Times New Roman")
plot(FIGDATA$Year,FIGDATA$spawn_dep+FIGDATA$tcb, pch=21,col="black", bg="blue",bty='L',ylab="Biomass (tons)", xlab="Year",cex.axis=1,cex.lab=1,
     ylim=c(0,MAXy), xaxt="n", family="A")
lines(FIGDATA$Year,FIGDATA$spawn_dep+FIGDATA$tcb,lty=2,lwd=1, col="black")
lines(FIGDATA$Year,FIGDATA$tot_mat_B_tons,lwd=3)
points(FIGDATA$Year,FIGDATA$spawn_dep+FIGDATA$tcb,pch=21, col="black",bg="blue", cex=1.3)
points(FIGDATA$Year,FIGDATA$tot_mat_B_tons,pch=17, col="red", cex=1.3)
points(max(FIGDATA$Year)+1,A$for_mat_B_st,pch=8,col="red",cex=1.3)#placeholder
textxy(max(FIGDATA$Year)+1, A$for_mat_B_st,round2(A$for_mat_B_st), pos=3, cex=1)
lines(FIGDATA$Year,FIGDATA$yminusone,lty=1,lwd=2, col="purple")
lines(FIGDATA$Year,FIGDATA$yminustwo,lty=1,lwd=2, col="green")
lines(FIGDATA$Year,FIGDATA$yminusthree,lty=1,lwd=2, col="orange")
#may need to change x-axis based on years used as past models
points(max(FIGDATA$Year),A$yminusoneFOR,pch=8,col="purple",cex=1.3)
points(max(FIGDATA$Year)-1,A$yminustwoFOR,pch=8,col="green",cex=1.3)
points(max(FIGDATA$Year)-2,A$yminusthreeFOR,pch=8,col="orange",cex=1.3)
textxy(max(FIGDATA$Year), A$yminusoneFOR,round2(A$yminusoneFOR), pos=1, cex=0.75)
textxy(max(FIGDATA$Year)-1, A$yminustwoFOR,round2(A$yminustwoFOR), pos=3, cex=0.75)
textxy(max(FIGDATA$Year)-2, A$yminusthreeFOR,round2(A$yminusthreeFOR), pos=1, cex=0.75)
lines(FIGDATA$Year,FIGDATA$Threshold,col="grey",lty=1, lwd=2)
axis(side=1,at=seq(min(FIGDATA$Year),max(FIGDATA$Year)+2,1),cex.axis=1, las=2)
legend("topleft",c("Survey estimated","Model estimated","2016 forecast","Threshold","2015 model estimated",
                   "2015 forecast","2014 model estimated","2014 forecast", "2013 model estimated","2013 forecast"), 
       lty=c(NA,NA,NA,1,1,NA,1,NA,1,NA),
       pch=c(16,17,8,NA, NA, 8, NA, 8, NA, 8),
       lwd=c(NA,NA,NA,2,2,NA,2,NA,2,NA), col=c("blue","red","red", "grey","purple",
                                               "purple", "green", "green", "orange", "orange"), cex=0.75, bty="n")
par(mar=c(5, 4, 2, 2) + 0.1)
dev.off()
C<-subset(FIGDATA, select=c(Year, tot_mat_B_tons, spawn_dep, tcb, yminusone, yminustwo, yminusthree))
Figure_3_data<-C
rm(MAXy,C)
#***********************************************************************************************************************
#GRAPHIC FOUR: 
#Figure 4: Residuals from model fits to survey egg deposition and Ricker spawner-recruit function.
detach(package:Hmisc, unload = TRUE) #messes with the barcharts in ggplot [fill="#E69F00"] 
attach(FIGDATA)
png(file='Figure 4.png', res=200, width=12, height=7, units ="in")  
g4a <- ggplot() +geom_bar(data=FIGDATA, mapping=aes(x=Year, y=res_tot_egg),stat='identity', position='dodge', fill="#E69F00") +
  ylab("Egg deposition residuals")
g4a<-g4a+scale_y_continuous(labels = fmt())+scale_x_continuous(breaks=seq(min(FIGDATA$Year),max(FIGDATA$Year)+1,1))+
  theme(axis.text.x = element_text(size=15,colour="black"),
        axis.title.x = element_text(size=15, colour="black"))+
  theme(axis.text.y = element_text(size=15,colour="black"),
        axis.title.y = element_text(size=15,colour="black"))+
theme(panel.border = element_rect(colour = "black"))+
  theme(legend.position="none")+
  theme(panel.grid.major = element_line(colour="white"))+
  theme(strip.text.x = element_text(size=14,face="bold"))
g4a<-g4a+theme_set(theme_bw(base_size=14,base_family=
                     'Times New Roman')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
cutoff <- data.frame(yintercept=0)
g4a<-g4a + geom_hline(aes(yintercept=yintercept), data=cutoff, show_guide=FALSE, colour="black", size=0.55)
g4a<-g4a+theme(axis.text.x=element_text(angle=-90))
B<-subset(FIGDATA, select=c(Year,res_SR)) 
B["Year"]<-B$Year+3
q<-max(B$Year)-2
B<-subset(B, B$Year<q)
B["Year"]<-as.numeric(B$Year)

g4b <- ggplot() +geom_bar(data=B, mapping=aes(x=Year, y=res_SR),stat='identity', position='dodge', fill="#009E73") +
  ylab("Spawner-recruit residuals")+xlab("Year")
g4b<-g4b+scale_y_continuous(labels = fmt())+theme_bw()+scale_x_continuous(breaks=seq(min(B$Year),max(B$Year+1),1))+
  theme(axis.text.x = element_text(size=14,colour="black"),
        axis.title.x = element_text(size=14, colour="black"))+
  theme(axis.text.y = element_text(size=14,colour="black"),
        axis.title.y = element_text(size=14,colour="black"))+
  theme(panel.border = element_rect(colour = "black"))+
  theme(legend.position="none")+
  theme(strip.text.x = element_text(size=14,face="bold"))
g4b<-g4b+theme_set(theme_bw(base_size=14,base_family=
                     'Times New Roman')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
cutoff <- data.frame(yintercept=0)
g4b<-g4b + geom_hline(aes(yintercept=yintercept), data=cutoff, show_guide=FALSE, colour="black", size=0.55)
g4b<-g4b+theme(axis.text.x=element_text(angle=-90))
library(grid)
g4<-arrange_ggplot2(g4a,g4b, ncol=1)
dev.off()
C<-subset(FIGDATA, select=c(Year, res_SR, res_tot_egg))
Figure_4_data<-C
rm(q,B,g4a,g4b,g4,C, cutoff)
ls(all.names = TRUE)
#***********************************************************************************************************************
#GRAPHIC FIVE: 
#Figure 5: Model estimates of age-3 recruit strength (numbers of age-3 mature and immature fish).
detach(package:Hmisc, unload = TRUE)

B<-subset(FIGDATA, select=c(Year, init_age_3)) 
maxY<-max(B$init_age_3,na.rm=TRUE)*1.5

g5 <- ggplot() +geom_bar(data=B, mapping=aes(x=Year, y=init_age_3),stat='identity', position='dodge', fill="#56B4E9") +
  ylab("Number of age-3 recruits (millions)")
g5<-g5 + theme(legend.position="none") +scale_x_continuous(breaks=seq(min(B$Year),max(B$Year+1),1))+coord_cartesian(ylim=c(0,maxY))
theme(axis.text.x = element_text(size=14, colour="black"),
        axis.title.x = element_text(size=14, colour="black"))+theme_bw()+
  theme(axis.text.y = element_text(size=14,colour="black"),
        axis.title.y = element_text(size=14,colour="black"))+
  theme(panel.background = element_rect(colour="white"))+
  theme(legend.position="none")+
  theme(panel.border = element_rect(colour = "black"))
g5<-g5+theme_set(theme_bw(base_size=14,base_family=
                            'Times New Roman')+
                   theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()))
g5<-g5+theme(axis.text.x=element_text(angle=-90))

png(file='Figure 5.png', res=200, width=10, height=6, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g5,vp=vplayout(1,1:1)) 
dev.off()
C<-subset(FIGDATA, select=c(Year, init_age_3))
Figure_5_data<-C
rm(B,g5,C, maxY)
ls(all.names = TRUE)
#***********************************************************************************************************************
#GRAPHIC SIX: 
#Figure 6: Spawning population biomass (blue bars;top figure), spawning population abundance (blue bars;middle figure), 
#population abundance (immature and spawning abundance) (blue bars;bottom figure), and commercial fishery 
#harvest (yellow bars) over time. The combination of the blue and yellow bars (total height of each bar) is the 
#mature biomass, mature population abundance, or total population abundance.
detach(package:Hmisc, unload = TRUE)
png(file='Figure 6.png', res=200, width=13, height=10, units ="in") 

B<-subset(FIGDATA, select=c(Year, tot_mat_B_tons, tot_sp_B_tons)) 
B["Catch"]<-B$tot_mat_B_tons-B$tot_sp_B_tons
B["postfishery"]<-B$tot_sp_B_tons

maxY<-max(B$postfishery,na.rm=TRUE)*1.5
B<-subset(B, select=c(Year,Catch,postfishery)) 
B<- melt(B, id=c("Year"), na.rm=TRUE)
B<- B[order(B$Year, -B$value) , ]
g6a <- ggplot() +geom_bar(data=B, mapping=aes(x=Year, y=value, fill=variable),stat="identity")+
  ylab("Biomass (tons)")+xlab("Year")+scale_fill_manual(values=c("#E69F00", "#0072B2"))
g6a<-g6a+theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())+  
                  theme(axis.text.x = element_text(size=12,colour="black", family="Times New Roman"),
                  axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
                  theme(axis.text.y = element_text(size=12,colour="black",family="Times New Roman"),
                         axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
                  theme(panel.background = element_rect(colour="white"))+
                  theme(panel.border = element_rect(colour = "black"))
g6a<-g6a+coord_cartesian(ylim=c(0,maxY))+scale_x_continuous(breaks=seq(min(B$Year),max(B$Year+1),1))+                  
     theme(legend.position="none")
g6a<-g6a+theme(axis.text.x=element_text(angle=-90))
rm(B,maxY)
ls(all.names = TRUE)

B<-subset(FIGDATA, select=c(Year, N, tot_post_N)) 
B["Catch"]<-B$N-B$tot_post_N
B["postfishery"]<-B$tot_post_N
maxY<-max(B$postfishery,na.rm=TRUE)*1.5
B<-subset(B, select=c(Year,Catch,postfishery)) 
B <- melt(B, id=c("Year"), na.rm=TRUE)
B<- B[order(B$Year, -B$value) , ]
g6b<-ggplot()+geom_bar(data=B, mapping=aes(x = Year,y=value, fill = variable),stat='identity') + ylab("Total abundance (millions)")+xlab("Year")+
scale_fill_manual(values=c("#E69F00", "#0072B2"))
g6b<-g6b+theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+  
  theme(axis.text.x = element_text(size=12,colour="black", family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))+
  theme(strip.text.x = element_text(size=14,face="bold", family="Times New Roman"))
g6b<-g6b+coord_cartesian(ylim=c(0,maxY))+scale_x_continuous(breaks=seq(min(B$Year),max(B$Year+1),1))+                  
  theme(legend.position="none")
g6b<-g6b+theme(axis.text.x=element_text(angle=-90))
rm(B,maxY)
ls(all.names = TRUE)

B<-subset(FIGDATA, select=c(Year, tot_sp_N, tot_mat_N)) 
B["Catch"]<-B$tot_mat_N-B$tot_sp_N
B["postfishery"]<-B$tot_sp_N
maxY<-max(B$postfishery,na.rm=TRUE)*1.5
B<-subset(B, select=c(Year,Catch,postfishery)) 
B<- melt(B, id=c("Year"), na.rm=TRUE)
B<- B[order(B$Year, -B$value) , ]
g6c<-ggplot()+geom_bar(data=B, mapping=aes(x = Year,y=value, fill = variable),stat='identity') + ylab("Abundance (millions)")+xlab("Year")+
  scale_fill_manual(values=c("#E69F00", "#0072B2"))
g6c<-g6c+theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+  
  theme(axis.text.x = element_text(size=12,colour="black", family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(panel.background = element_rect(colour="white"))+
  theme(panel.border = element_rect(colour = "black"))+
  theme(strip.text.x = element_text(size=14,face="bold", family="Times New Roman"))
g6c<-g6c+coord_cartesian(ylim=c(0,maxY))+scale_x_continuous(breaks=seq(min(B$Year),max(B$Year+1),1))+                  
  theme(legend.position="none")
g6c<-g6c+theme(axis.text.x=element_text(angle=-90))
library(grid)
g6<-arrange_ggplot2(g6a,g6c,g6b, ncol=1)
dev.off()
C<-subset(FIGDATA, select=c(Year, tot_mat_B_tons, tot_sp_B_tons, N, tot_post_N, tot_sp_N, tot_mat_N))
Figure_6_data<-C
rm(g6,g6a,g6b,B,g6c,maxY,C)
ls(all.names = TRUE)
#***********************************************************************************************************************
#GRAPHIC SEVEN: 
#Figure 7: Model estimates of gear selectivity at age (a), maturity at age (b), and survival (c) by year.
#Figure 7a: Model estimates of gear selectivity at age.
detach(package:Hmisc, unload = TRUE)
y<-as.data.frame(A$Scaled_Gear_Selectivity)
B<-subset(FIGDATA, select=c(Year)) #predicted
C <- cbind(y,B)
C <- merge(FIGDATA,C,by=c("Year"), all=TRUE)
C["Age3"]<-C$V1
C["Age4"]<-C$V2
C["Age5"]<-C$V3
C["Age6"]<-C$V4
C["Age7"]<-C$V5
C["Age8"]<-C$V6
C<-subset(C, select=c(Year, Age3, Age4, Age5, Age6, Age7, Age8))
C<- melt(C, id=c("Year"), na.rm=TRUE)
C["Age"] <- ifelse(C$variable=="Age2","2", ifelse (C$variable=="Age3","3",
                                                   ifelse (C$variable=="Age4","4",
                                                           ifelse (C$variable=="Age5","5",
                                                                   ifelse (C$variable=="Age6","6",
                                                                           ifelse (C$variable=="Age7","7","8+"))))))
C["GS"]<-C$value
C<-subset(C, select=c(Year, Age, GS))
g7a<-ggplot(data=C,aes(x=Age, y=GS,fill=Age))+facet_wrap(~Year,ncol=7,as.table=TRUE)+geom_bar(stat="identity")+
  geom_line(data=C,aes(x=Age,y=GS,group=Year),show_guide=FALSE,size=1,colour="black")+ theme_bw()+
   xlab ("Age")+ ylab("Proportion-at-age")+
   theme(text=element_text(family="Times New Roman", face="bold", size=12))
g7a<-g7a+coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(plot.title=element_text(size=rel(1.5),colour="black",vjust =1))+
  theme(strip.text.x = element_text(size=14,face="bold", family="Times New Roman"))+
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())


png(file='Figure 7a.png', res=200, width=9, height=6, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g7a,vp=vplayout(1,1:1)) 
Figure_7a_data<-C
dev.off()
rm(B,C,y,g7a)   
ls(all.names = TRUE)

#Figure 7b: Model estimates of maturity-at-age.
y<-as.data.frame(A$Maturity)
B<-subset(FIGDATA, select=c(Year)) #predicted
C <- cbind(y,B)
C <- merge(FIGDATA,C,by=c("Year"), all=TRUE)
C["Age3"]<-C$V1
C["Age4"]<-C$V2
C["Age5"]<-C$V3
C["Age6"]<-C$V4
C["Age7"]<-C$V5
C["Age8"]<-C$V6
C<-subset(C, select=c(Year, Age3, Age4, Age5, Age6, Age7, Age8))
C<- melt(C, id=c("Year"), na.rm=TRUE)
C["Age"] <- ifelse(C$variable=="Age2","2", ifelse (C$variable=="Age3","3",
                                                   ifelse (C$variable=="Age4","4",
                                                           ifelse (C$variable=="Age5","5",
                                                                   ifelse (C$variable=="Age6","6",
                                                                           ifelse (C$variable=="Age7","7","8+"))))))
C["MAT"]<-C$value
C<-subset(C, select=c(Year, Age, MAT))
g7b<-ggplot(data=C,aes(x=Age, y=MAT,fill=Age))+facet_wrap(~Year,ncol=7,as.table=TRUE)+geom_bar(stat="identity")+
  geom_line(data=C,aes(x=Age,y=MAT,group=Year),show_guide=FALSE,size=1,colour="black")+ theme_bw()+
 xlab ("Age")+ ylab("Proportion-at-age")+theme(text=element_text(family="Times New Roman", face="bold", size=12))
g7b<-g7b+coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(size=12,colour="black", family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black", family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black", family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black", family="Times New Roman"))+
  theme(plot.title=element_text(size=rel(1.5),colour="black",vjust =1))+
  theme(strip.text.x = element_text(size=14,face="bold", family="Times New Roman"))+
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
png(file='Figure 7b.png', res=200, width=9, height=6, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g7b,vp=vplayout(1,1:1)) 
dev.off()
Figure_7b_data<-C
rm(y,B,C,g7b)   
ls(all.names = TRUE)

#Figure 7c: Model estimates of Survival.
y<-as.data.frame(A$Sur)
B<-subset(FIGDATA, select=c(Year)) #predicted
C <- cbind(y,B)
C <- merge(FIGDATA,C,by=c("Year"), all=TRUE)
C["Age3"]<-C$V1
C["Age4"]<-C$V2
C["Age5"]<-C$V3
C["Age6"]<-C$V4
C["Age7"]<-C$V5
C["Age8"]<-C$V6
C<-subset(C, select=c(Year, Age3, Age4, Age5, Age6, Age7, Age8))
C<- melt(C, id=c("Year"), na.rm=TRUE)
C["Age"] <- ifelse(C$variable=="Age2","2", ifelse (C$variable=="Age3","3",
                                                   ifelse (C$variable=="Age4","4",
                                                           ifelse (C$variable=="Age5","5",
                                                                   ifelse (C$variable=="Age6","6",
                                                                           ifelse (C$variable=="Age7","7","8+"))))))
C["S"]<-C$value
C<-subset(C, select=c(Year, Age, S))
g7c<-ggplot(data=C,aes(x=Age, y=S,fill=Age))+facet_wrap(~Year,ncol=7,as.table=TRUE)+geom_bar(stat="identity")+
  geom_line(data=C,aes(x=Age,y=S,group=Year),show_guide=FALSE,size=1,colour="black")+ theme_bw()+
xlab ("Age")+ylab("Proportion-at-age")+theme(text=element_text(family="Times New Roman", face="bold", size=12))
g7c<-g7c+coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(plot.title=element_text(size=rel(1.5),colour="black",vjust =1))+
  theme(strip.text.x = element_text(size=14,face="bold",family="Times New Roman"))+ 
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
png(file='Figure 7c.png', res=200, width=9, height=6, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g7c,vp=vplayout(1,1:1)) 
dev.off()
Figure_7c_data<-C
rm(y,B,C,g7c)  
ls(all.names = TRUE)
#***********************************************************************************************************************
#GRAPHIC EIGHT: 
#Figure 8:Observed seine or gillnet (red line with square points) and model-estimated (bar) catch-age composition.
detach(package:Hmisc, unload = TRUE)
D<-subset(FIGDATA, select=c(Year,sel_naa_prop3, sel_naa_prop4, sel_naa_prop5, sel_naa_prop6, sel_naa_prop7, sel_naa_prop8)) #estimated
D["Age3"]<-D$sel_naa_prop3
D["Age4"]<-D$sel_naa_prop4
D["Age5"]<-D$sel_naa_prop5
D["Age6"]<-D$sel_naa_prop6
D["Age7"]<-D$sel_naa_prop7
D["Age8"]<-D$sel_naa_prop8
D["SUM"]<-sum(D$Age3, D$Age4, D$Age5, D$Age6, D$Age7, D$Age8)
D<-subset(D, select=c(Year,Age3, Age4, Age5, Age6, Age7, Age8))
D<- melt(D, id=c("Year"), na.rm=TRUE)
D["Age"] <- ifelse(D$variable=="Age2","2", ifelse (D$variable=="Age3","3",
                                                   ifelse (D$variable=="Age4","4",
                                                           ifelse (D$variable=="Age5","5",
                                                                   ifelse (D$variable=="Age6","6",
                                                                           ifelse (D$variable=="Age7","7","8+"))))))

D["EST"]<-D$value
D<-subset(D, select=c(Year,Age, EST))
B<-subset(FIGDATA, select=c(Year,obs_c_comp3, obs_c_comp4, obs_c_comp5, obs_c_comp6, obs_c_comp7, obs_c_comp8))#observed
B["Age3"]<-B$obs_c_comp3
B["Age4"]<-B$obs_c_comp4
B["Age5"]<-B$obs_c_comp5
B["Age6"]<-B$obs_c_comp6
B["Age7"]<-B$obs_c_comp7
B["Age8"]<-B$obs_c_comp8
B["SUM"]<-sum(B$Age3, B$Age4, B$Age5, B$Age6, B$Age7, B$Age8)
B<-subset(B, select=c(Year,Age3, Age4, Age5, Age6, Age7, Age8))
B<- melt(B, id=c("Year"), na.rm=TRUE)
B["Age"] <- ifelse(B$variable=="Age2","2", ifelse (B$variable=="Age3","3",
                                                   ifelse (B$variable=="Age4","4",
                                                           ifelse (B$variable=="Age5","5",
                                                                   ifelse (B$variable=="Age6","6",
                                                                           ifelse (B$variable=="Age7","7","8+"))))))

B["OBS"]<-B$value
B<-subset(B, select=c(Year,Age, OBS))
C <- merge(D,B,by=c("Year", "Age"), all=TRUE)

#Sitka facet_wrap(~Year,nrow=7, ncol=6)
C <- rbind(C,data.frame(Year=2016,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2017,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2018,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2019,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2020,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2021,Age=3,EST=NA, OBS=NA))
C$Year<-factor(C$Year,levels=c(
"1980", "1987","1994","2001", "2008","2015",
"1981","1988","1995","2002", "2009","2016",
"1982","1989","1996","2003", "2010", "2017",
"1983","1990","1997","2004","2011", "2018",
"1984","1991","1998","2005","2012","2019",
"1985","1992","1999","2006", "2013","2020",
"1986","1993","2000","2007","2014","2021"))

#ncol & nrow needs to be changed based on stock*****************
g8<-ggplot(data=C,aes(x=Age, y=EST,fill=Age))+facet_wrap(~Year,nrow=7, ncol=6)+geom_bar(stat="identity")+
  geom_line(data=C,aes(x=Age,y=OBS,group=Year),show_guide=FALSE,size=1,colour="red")+
  geom_point(data = C,aes(x=Age,y=OBS,group=Year),show_guide=FALSE,size=2,shape=22,fill="red")+ theme_bw()+
 xlab ("Age")+ylab("Proportion-at-age")+theme(text=element_text(family="Times New Roman", face="bold", size=12))
g8<-g8+coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(plot.title=element_text(size=rel(1.5),colour="black",vjust =1,family="Times New Roman"))+
  theme(strip.text.x = element_text(size=14,face="bold"))+
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
png(file='Figure 8.png', res=200, width=9, height=9, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g8,vp=vplayout(1,1:1)) 
dev.off()
Figure_8_data<-C
rm(B,C,g8,D)
ls(all.names = TRUE)
#*************************************************************************************************************************
#GRAPHIC NINE: 
#Figure 9:Observed cast net (red line with square points) and model-estimated (bar) spawning-age composition.
detach(package:Hmisc, unload = TRUE)
D<-subset(FIGDATA, select=c(Year,est_sp_comp3, est_sp_comp4, est_sp_comp5, est_sp_comp6, est_sp_comp7, est_sp_comp8)) #estimated
D["Age3"]<-D$est_sp_comp3
D["Age4"]<-D$est_sp_comp4
D["Age5"]<-D$est_sp_comp5
D["Age6"]<-D$est_sp_comp6
D["Age7"]<-D$est_sp_comp7
D["Age8"]<-D$est_sp_comp8
D["SUM"]<-sum(D$Age3, D$Age4, D$Age5, D$Age6, D$Age7, D$Age8)
D<-subset(D, select=c(Year,Age3, Age4, Age5, Age6, Age7, Age8))
D<- melt(D, id=c("Year"), na.rm=TRUE)
D["Age"] <- ifelse(D$variable=="Age2","2", ifelse (D$variable=="Age3","3",
                                                    ifelse (D$variable=="Age4","4",
                                                            ifelse (D$variable=="Age5","5",
                                                                    ifelse (D$variable=="Age6","6",
                                                                            ifelse (D$variable=="Age7","7","8+"))))))

D["EST"]<-D$value
D<-subset(D, select=c(Year,Age, EST))
B<-subset(FIGDATA, select=c(Year,obs_sp_comp3, obs_sp_comp4, obs_sp_comp5, obs_sp_comp6, obs_sp_comp7, obs_sp_comp8)) #observed
B["Age3"]<-B$obs_sp_comp3
B["Age4"]<-B$obs_sp_comp4
B["Age5"]<-B$obs_sp_comp5
B["Age6"]<-B$obs_sp_comp6
B["Age7"]<-B$obs_sp_comp7
B["Age8"]<-B$obs_sp_comp8
B["SUM"]<-sum(B$Age3, B$Age4, B$Age5, B$Age6, B$Age7, B$Age8)
B<-subset(B, select=c(Year,Age3, Age4, Age5, Age6, Age7, Age8))
B<- melt(B, id=c("Year"), na.rm=TRUE)
B["Age"] <- ifelse(B$variable=="Age2","2", ifelse (B$variable=="Age3","3",
                                                    ifelse (B$variable=="Age4","4",
                                                            ifelse (B$variable=="Age5","5",
                                                                    ifelse (B$variable=="Age6","6",
                                                                            ifelse (B$variable=="Age7","7","8+"))))))

B["OBS"]<-B$value
B<-subset(B, select=c(Year,Age, OBS))
C <- merge(D,B,by=c("Year", "Age"), all=TRUE)
#only run for code selected stock*****************
#Sitka facet_wrap(~Year,nrow=7, ncol=6)
C <- rbind(C,data.frame(Year=2016,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2017,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2018,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2019,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2020,Age=3,EST=NA, OBS=NA))
C <- rbind(C,data.frame(Year=2021,Age=3,EST=NA, OBS=NA))
C$Year<-factor(C$Year,levels=c(
  "1980", "1987","1994","2001", "2008","2015",
  "1981","1988","1995","2002", "2009","2016",
  "1982","1989","1996","2003", "2010", "2017",
  "1983","1990","1997","2004","2011", "2018",
  "1984","1991","1998","2005","2012","2019",
  "1985","1992","1999","2006", "2013","2020",
  "1986","1993","2000","2007","2014","2021"))

#ncol & nrow needs to be changed based on stock*****************
g9<-ggplot(data=C,aes(x=Age, y=EST,fill=Age))+facet_wrap(~Year,nrow=7, ncol=6)+geom_bar(stat="identity")+
  geom_line(data=C,aes(x=Age,y=OBS,group=Year),show_guide=FALSE,size=1,colour="red")+
  geom_point(data = C,aes(x=Age,y=OBS,group=Year),show_guide=FALSE,size=2,shape=22,fill="red")+ theme_bw()+
xlab("Age")+ylab("Proportion-at-age")+theme(text=element_text(family="Times New Roman", face="bold", size=12))
g9<-g9+coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(size=12,colour="black", family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black", family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black", family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black", family="Times New Roman"))+
  theme(plot.title=element_text(size=rel(1.5),colour="black",vjust =1, family="Times New Roman"))+
  theme(strip.text.x = element_text(size=14,face="bold", family="Times New Roman"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
png(file='Figure 9.png', res=200, width=9, height=9, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g9,vp=vplayout(1,1:1)) 
Figure_9_data<-C
dev.off()
rm(D,B,C,g9)
#*************************************************************************************************************************
#GRAPHIC TEN: 
#Figure 10: Spawning biomass (tons) versus age-3 abundance (millions of mature and immature fish) (blue circles) 
#with Ricker-estimated age-3 abundance (red triangles), from 1991-2014. 
D<-subset(FIGDATA, select=c(Year, tot_sp_B_tons, SR)) 
D["Year"]<-D$Year+3
D<-subset(D, D$Year<max(D$Year)-2)
B<-subset(FIGDATA, select=c(Year, init_age_3)) 
B<-subset(B, B$Year>min(B$Year)+2)
C <- merge(D, B, by=c("Year"), all=TRUE) #Merges datasets
MAXy<-max(C$init_age_3,na.rm=TRUE)*1.5 # set y max. limit

png(file='Figure 10.png', res=200, width=7, height=5, units ="in")  
windowsFonts(A = windowsFont("Times New Roman"))
op <- par(family = "Times New Roman")
getOption("scipen")
opt <- options("scipen" = 20)
getOption("scipen")
plot(C$tot_sp_B_tons,C$init_age_3,pch=21,col="black",bg="blue",cex=1,ylab="Number of recruits (millions)",
     xlab="Spawning biomass (tons)",cex.axis=1,cex.lab=1.2, ylim=c(0, MAXy), family="A")
points(C$tot_sp_B_tons,C$SR,pch=24,col="black", bg="red", cex=1) #predicted
textxy(C$tot_sp_B_tons, C$init_age_3, C$Year, pos=3, cex=0.7)
options(opt)
legend("topright", c("ASA Age-3 estimated abundance","Ricker estimated"),pch=c(16,17),col=c("blue", "red"),
       ,cex=1, bty="n")
dev.off()
Figure_10_data<-C
rm (D,B,C, MAXy)
#***********************************************************************************************************************
#GRAPHIC ELEVEN:
#Figure 11: Projected mature biomass at age (tons) for forecast year.
detach(package:Hmisc, unload = TRUE)
MAXy<-max(FIGDATAAGE$for_mat_baa_tons,na.rm=TRUE)*1.5 # set y max. limit
FIGDATAAGE["for_mat_baa_tons"]<-as.numeric(FIGDATAAGE$for_mat_baa_tons)
#Check to make sure that the numbers on top of bar sum to to forecasted biomass
#If they don't enter the % manually below
#FIGDATAAGE["for_mat_baa_tons"]<-c(1578,3480,157,1589,882,1566)
g11<-ggplot(FIGDATAAGE, aes(x=Age2, y=for_mat_baa_tons)) +  ylab("Forecasted Biomass (tons)")+xlab("Age")+
  theme(text=element_text(family="Times New Roman", size=12))+
  geom_bar(stat="identity",position="dodge" ,fill="#E69F00")

g11<-g11+geom_text(data=FIGDATAAGE,aes(label=round2(for_mat_baa_tons), cex=1),
                   vjust=-1)
g11<-g11+coord_cartesian(ylim=c(0,MAXy))
g11<-g11 + theme(legend.position="none")+
  theme(axis.text.x = element_text(size=14,colour="black",family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=14,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(strip.text.x = element_text(size=14,family="Times New Roman"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
png(file='Figure 11.png', res=200, width=6, height=4, units ="in")   
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g11,vp=vplayout(1,1:1)) 
dev.off()
C<-subset(FIGDATAAGE, select=c(Age2, for_mat_baa_tons))
Figure_11_data<-C
rm(MAXy,g11)
ls(all.names = TRUE)
#***********************************************************************************************************************
#GRAPHIC TWELVE: *ADD % to top of graph
#Figure 12: Projected percentage of mature numbers at age for forecast year. 
#*Check to make sure that the numbers on top of bar sum to 100!!!!
detach(package:Hmisc, unload = TRUE)
D<-subset(FIGDATAAGE, select=c(for_mat_prop)) 
D["Percentage"]<-round(D$for_mat_prop,2)
D<-subset(D, select=c(Percentage)) 
x<-(colSums(D))
D["Percentage1"]<-round(D$Percentage/x,2)
D["Percentage2"]<-D$Percentage1
D<-subset(D, select=c(Percentage2))
Check<-sum(D$Percentage2)#Check to make sure that the numbers on top of bar sum to 100!!!!
#If they don't enter the % manually below
#D["Percentage2"]<-c(19,7,32,15,16,11)
Check<-sum(D$Percentage2)
Check
C<-cbind(FIGDATAAGE, D)
C<-cbind(FIGDATAAGE, D)#Check to make sure that the numbers on top of bar sum to 100!!!!
C$Age <- factor(C$Age)
g12<-ggplot(C, aes(x=Age2, y=Percentage2*100, fill=Age)) + ylab("Forecasted Percentage")+xlab("Age")+theme(text=element_text(family="Times New Roman", size=12))+
  geom_bar(stat="identity",position="dodge")
g12<-g12+scale_fill_manual(values=c("#FF0000","#3A5FCD", "#3A5FCD", "#3A5FCD","#3A5FCD", "#3A5FCD"))
g12<-g12+coord_cartesian(ylim=c(0,100))+geom_text(data=C,aes(label=paste(C$Percentage2*100, "%", sep=""),cex=1),vjust=-1, family="Times New Roman") + #add prop. on top of bars
  theme(axis.text.x = element_text(size=14,colour="black",family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=14,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(legend.position="none")+
  theme(strip.text.x = element_text(size=14,face="bold",family="Times New Roman"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

png(file='Figure 12.png', res=200, width=6, height=4, units ="in")   
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g12,vp=vplayout(1,1:1)) 
dev.off()
C<-subset(FIGDATAAGE, select=c(Age2, for_mat_prop))
Figure_12_data<-C
rm(g12, D, x, C, Check)
ls(all.names = TRUE)
#***********************************************************************************************************************
#GRAPHIC THIRTEEN: 
#Figure 13: Forecasted weight-at-age. 
detach(package:Hmisc, unload = TRUE)
MAXy<-max(FIGDATAAGE$fw_a_a,na.rm=TRUE)*1.5 # set y max. limit
FIGDATAAGE["fw_a_a"]<-round(FIGDATAAGE$fw_a_a,1)
g13<-ggplot(FIGDATAAGE, aes(x=Age2, y=fw_a_a)) +  ylab("Forecasted Weight (g)")+xlab("Age")+theme(text=element_text(family="Times New Roman", size=12))+
  geom_bar(stat="identity",position="dodge" ,fill="#56B4E9") 
g13<-g13+geom_text(data=FIGDATAAGE,aes(label=paste(FIGDATAAGE$fw_a_a,"g", sep=""), cex=1),vjust=-1) #add prop. on top of bars
g13<-g13+coord_cartesian(ylim=c(0,MAXy))
g13<-g13 + theme(legend.position="none")+
theme(axis.text.x = element_text(size=14,colour="black",family="Times New Roman"),
      axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=14,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(strip.text.x = element_text(size=14,face="bold",family="Times New Roman"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
png(file='Figure 13.png', res=200, width=6, height=4, units ="in")   
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g13,vp=vplayout(1,1:1)) 
dev.off()
C<-subset(FIGDATAAGE, select=c(Age2, fw_a_a))
Figure_13_data<-C
rm(MAXy,g13,C)
ls(all.names = TRUE)
#***********************************************************************************************************************
#GRAPHIC FOURTEEN: 
#Figure 14: Stacked bar graph of catch (orange), spawning biomass (green), GHL (blue), and the spawning biomass forecast (pink). The harvest (or GHL) plus the spawning biomass equals the mature biomass. 
detach(package:Hmisc, unload = TRUE)
#The pre-fishery-postfishery=catch
D<-subset(FIGDATA, select=c(Year, tot_mat_B_tons, tot_sp_B_tons)) 
D["Catch"]<-D$tot_mat_B_tons-D$tot_sp_B_tons
D["Spawning_biomass"]<-D$tot_sp_B_tons 
D["GHL"]<-0
D["Spawning_biomass_forecast"]<-0
D<-subset(D, select=c(Year,Catch,Spawning_biomass, GHL, Spawning_biomass_forecast))
maxY<-max(D$Spawning_biomass,na.rm=TRUE)*1.2

for_mat_B_st<-A$for_mat_B_st
GHL<-A$GHL
Year<-c(max(D$Year)+1)
Catch<-0
Spawning_biomass<-0 
Spawning_biomass_forecast<-for_mat_B_st-GHL 
B<- data.frame(Year, Catch, Spawning_biomass,GHL, Spawning_biomass_forecast) 
C <- rbind(D,B)
C <- melt(C, id=c("Year"), na.rm=TRUE)
C$variable <- factor(C$variable, levels=c("Catch", "Spawning_biomass","GHL","Spawning_biomass_forecast"))

cbPalette <- c("#E69F00","#0072B2","#009E73","#CC79A7")
C["variable"]<-ifelse(C$variable=="Catch","Catch", ifelse (C$variable=="GHL","GHL",
                          ifelse (C$variable=="Spawning_biomass","Spawning biomass","Spawning biomass forecast")))

  
  
C<- C[order(C$Year, -C$value) , ]
g14 <- ggplot() +geom_bar(data=C, mapping=aes(x=Year, y=value, fill=variable),stat='identity')+ scale_fill_manual(values=cbPalette)+
ylab("Biomass (tons)")+xlab("Year")+theme(text=element_text(family="Times New Roman", size=12))
g14<-g14 +coord_cartesian(ylim=c(0, maxY)) + 
  scale_y_continuous(breaks=seq(0, maxY, 5000))+theme_bw()+scale_x_continuous(breaks=seq(min(D$Year),max(D$Year+2),1))+
  theme(axis.text.x = element_text(size=12,colour="black",family="Times New Roman"),
      axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(panel.border = element_rect(colour = "black"))+
  theme(legend.title=element_blank())+
  theme(strip.text.x = element_text(size=14,face="bold",family="Times New Roman"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position=c(0.2,0.9))
g14<-g14+theme(axis.text.x=element_text(angle=-90))

x<-subset(FIGDATA, select=c(Year, Threshold))#estimated mature biomass
Y<-c(max(x$Year)+1)
B<- data.frame(Y)  
B["Year"]<-B$Y
B["Threshold"]<-A$Thresh
B<-subset(B, select=c(Year, Threshold))
B<-rbind(x,B) #to create x axis with a greater x value than max(Year)
addlinetoplot <- function(dataset, varx, vary) { 
  list(
    geom_line(data=dataset, aes_string(x=varx, y=vary)), 
    geom_point(data=dataset, aes_string(x=varx, y=vary))
  )}

g14<-g14 + geom_line(data=B,aes(x=Year, y=Threshold), show_guide=FALSE, colour="grey", size=2)
png(file='Figure 14.png', res=200, width=12, height=8, units ="in") 
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g14,vp=vplayout(1,1:1)) 
dev.off()
Figure_14_data<-C
rm(B,C, D, g14,maxY,Spawning_biomass_forecast,Spawning_biomass,Year,GHL,Catch)
ls(all.names = TRUE)
#************************************************************************************************************************
#GRAPHIC FIFTEEN: 
#Figure 15: Spawning Age Compostion residuals.
detach(package:Hmisc, unload = TRUE)
x<-as.data.frame(A$res_sp_comp)
x["Age3"]<-x$V1
x["Age4"]<-x$V2
x["Age5"]<-x$V3
x["Age6"]<-x$V4
x["Age7"]<-x$V5
x["Age8"]<-x$V6
y<-subset(FIGDATA, select=c(Year)) 
C <-cbind(x,y)
C<-subset(C, select=c(Year, Age3, Age4, Age5,Age6,Age7, Age8)) 
C<- melt(C, id=c("Year"), na.rm=TRUE)
C["Age"] <- ifelse(C$variable=="Age2","2", ifelse (C$variable=="Age3","3",
                                                   ifelse (C$variable=="Age4","4",
                                                           ifelse (C$variable=="Age5","5",
                                                                   ifelse (C$variable=="Age6","6",
                                                                           ifelse (C$variable=="Age7","7","8+"))))))

C["FREQ"]<-C$value
C["Variable"]<-"Estimated"
C<-subset(C, select=c(Year,Age, FREQ, Variable))

g15<-ggplot(data=C,aes(x=Age, y=FREQ, fill=Variable))+facet_wrap(~Year,nrow=6, ncol=6)+geom_bar(stat="identity", position="dodge")+
 theme_bw()+xlab ("Age")+ylab("Proportion-at-age")+
  scale_fill_manual(values=c("#E69F00", "#0072B2"))+theme(legend.title=element_blank())+theme(text=element_text(family="Times New Roman", size=12))
g15<-g15+ theme(axis.text.x = element_text(size=12,colour="black", family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black", family="Times New Roman"))+
  theme(axis.text.y = element_text(size=12,colour="black", family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black" ,family="Times New Roman"))+
  theme(plot.title=element_text(size=rel(1.5),colour="black",vjust =1,family="Times New Roman" ))+
  theme(strip.text.x = element_text(size=14,face="bold", family="Times New Roman"))+theme(legend.position="none")+
 theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

png(file='Figure 15.png', res=200, width=9, height=9, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g15,vp=vplayout(1,1:1)) 
dev.off()
Figure_15_data<-C
rm(x,y,C,g15)
#************************************************************************************************************************
#GRAPHIC SIXTEEN: 
#Figure 16: Spawning age composition figure for model analysis only.
my.col <- colorRampPalette(brewer.pal(9,"RdYlGn"))
# THESE COMMANDS CAN BE RUN AND THE COLOR SCHEME MODIFIED; IF SO, THEN I'D 
# RECOMMEND EXPLICITLY INCLUDING THE CODE AS WRITTEN SO THE MODS DON'T 
# HAVE TO BE REDONE
#plot.table.helper.color <- edit(plot.table.helper.color)
#plot.table.helper.colorbar <- edit(plot.table.helper.colorbar)
x<-as.data.frame(A$res_sp_comp, round=4)
x<-round(x,4)
options(scipen=999) 
x$Year<-FIGDATA$Year

temp<-as.matrix(x)

colnames(temp)<-c(paste('Age',3:7,sep=' '), "Age 8+","Year")
temp<-subset(temp,select=c("Year", "Age 3", "Age 4", "Age 5","Age 6","Age 7", "Age 8+"))

# plot temp with colorbar, display Correlation in (top, left) cell
png(file='Figure 16.png', res=200, width=9, height=9, units ="in")  
plot.table(temp, smain='Residual', highlight = TRUE, colorbar = TRUE)

dev.off()
rm(x)
#************************************************************************************************************************
#GRAPHIC SEVENTEEN 
#Figure 17: Catch age composition figure for model analysis only.
my.col <- colorRampPalette(brewer.pal(9,"RdYlGn"))

# THESE COMMANDS CAN BE RUN AND THE COLOR SCHEME MODIFIED; IF SO, THEN I'D 
# RECOMMEND EXPLICITLY INCLUDING THE CODE AS WRITTEN SO THE MODS DON'T 
# HAVE TO BE REDONE
#plot.table.helper.color <- edit(plot.table.helper.color)
#plot.table.helper.colorbar <- edit(plot.table.helper.colorbar)
x<-as.data.frame(A$res_c_comp)
x<-round(x,4)
options(scipen=999) 
x$Year<-FIGDATA$Year

temp<-as.matrix(x)

colnames(temp)<-c(paste('Age',3:7,sep=' '), "Age 8+","Year")
temp<-subset(temp,select=c("Year", "Age 3", "Age 4", "Age 5","Age 6","Age 7", "Age 8+"))

# plot temp with colorbar, display Correlation in (top, left) cell
png(file='Figure 17.png', res=200, width=9, height=9, units ="in")  
plot.table(temp, smain='Residual', highlight = TRUE, colorbar = TRUE)

dev.off()
rm(x)
#GRAPHIC EIGHTa: 
#Figure 8a:Observed seine or gillnet (red line with square points) and model-estimated (bar) catch-age composition.
detach(package:Hmisc, unload = TRUE)
D<-subset(FIGDATA, select=c(Year,sel_naa_prop3, sel_naa_prop4, sel_naa_prop5, sel_naa_prop6, sel_naa_prop7, sel_naa_prop8)) #estimated
D["Age3"]<-D$sel_naa_prop3
D["Age4"]<-D$sel_naa_prop4
D["Age5"]<-D$sel_naa_prop5
D["Age6"]<-D$sel_naa_prop6
D["Age7"]<-D$sel_naa_prop7
D["Age8"]<-D$sel_naa_prop8
D["SUM"]<-sum(D$Age3, D$Age4, D$Age5, D$Age6, D$Age7, D$Age8)
D<-subset(D, select=c(Year,Age3, Age4, Age5, Age6, Age7, Age8))
D<- melt(D, id=c("Year"), na.rm=TRUE)
D["Age"] <- ifelse(D$variable=="Age2","2", ifelse (D$variable=="Age3","3",
                                                   ifelse (D$variable=="Age4","4",
                                                           ifelse (D$variable=="Age5","5",
                                                                   ifelse (D$variable=="Age6","6",
                                                                           ifelse (D$variable=="Age7","7","8+"))))))

D["EST"]<-D$value
D["Variable"]<-"Estimated"
D<-subset(D, select=c(Year,Age, EST, Variable))
B<-subset(FIGDATA, select=c(Year,obs_c_comp3, obs_c_comp4, obs_c_comp5, obs_c_comp6, obs_c_comp7, obs_c_comp8))#observed
B["Age3"]<-B$obs_c_comp3
B["Age4"]<-B$obs_c_comp4
B["Age5"]<-B$obs_c_comp5
B["Age6"]<-B$obs_c_comp6
B["Age7"]<-B$obs_c_comp7
B["Age8"]<-B$obs_c_comp8
B["SUM"]<-sum(B$Age3, B$Age4, B$Age5, B$Age6, B$Age7, B$Age8)
B<-subset(B, select=c(Year,Age3, Age4, Age5, Age6, Age7, Age8))
B<- melt(B, id=c("Year"), na.rm=TRUE)
B["Age"] <- ifelse(B$variable=="Age2","2", ifelse (B$variable=="Age3","3",
                                                   ifelse (B$variable=="Age4","4",
                                                           ifelse (B$variable=="Age5","5",
                                                                   ifelse (B$variable=="Age6","6",
                                                                           ifelse (B$variable=="Age7","7","8+"))))))

B["EST"]<-B$value
B["Variable"]<-"Observed"
B<-subset(B, select=c(Year,Age, EST, Variable))
C <- rbind(B,D)
C <- rbind(C,data.frame(Year=2016,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2017,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2018,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2019,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2020,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2021,Age=3,EST=NA, Variable="Estimated"))
C$Year<-factor(C$Year,levels=c(
  "1980", "1987","1994","2001", "2008","2015",
  "1981","1988","1995","2002", "2009","2016",
  "1982","1989","1996","2003", "2010", "2017",
  "1983","1990","1997","2004","2011", "2018",
  "1984","1991","1998","2005","2012","2019",
  "1985","1992","1999","2006", "2013","2020",
  "1986","1993","2000","2007","2014","2021"))
g8a<-ggplot(data=C,aes(x=Age, y=EST, fill=Variable))+facet_wrap(~Year,nrow=7, ncol=6)+geom_bar(stat="identity", position="dodge")+
  theme_bw()+xlab ("Age")+ylab("Proportion-at-age")+
  scale_fill_manual(values=c("#E69F00", "#0072B2"))+theme(legend.title=element_blank())+theme(text=element_text(family="Times New Roman", face="bold", size=12))
g8a<-g8a+coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(size=10,colour="black", family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black",family="Times New Roman"))+
  theme(axis.text.y = element_text(size=10,colour="black",family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black",family="Times New Roman"))+
  theme(plot.title=element_text(size=rel(1.5),colour="black",vjust =1,family="Times New Roman"))+
  theme(strip.text.x = element_text(size=14,face="bold",family="Times New Roman"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
png(file='Figure 8a.png', res=200, width=9, height=9, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g8a,vp=vplayout(1,1:1)) 
dev.off()
Figure_8a_data<-C
rm(B,g8a,D,C)
#***********************************************************************************************************************
#***********************************************************************************************************************
#GRAPHIC NINEa: 
#Figure 9a:Observed cast net (red line with square points) and model-estimated (bar) spawning-age composition.
detach(package:Hmisc, unload = TRUE)
D<-subset(FIGDATA, select=c(Year,est_sp_comp3, est_sp_comp4, est_sp_comp5, est_sp_comp6, est_sp_comp7, est_sp_comp8)) #estimated
D["Age3"]<-D$est_sp_comp3
D["Age4"]<-D$est_sp_comp4
D["Age5"]<-D$est_sp_comp5
D["Age6"]<-D$est_sp_comp6
D["Age7"]<-D$est_sp_comp7
D["Age8"]<-D$est_sp_comp8
D["SUM"]<-sum(D$Age3, D$Age4, D$Age5, D$Age6, D$Age7, D$Age8)
D<-subset(D, select=c(Year,Age3, Age4, Age5, Age6, Age7, Age8))
D<- melt(D, id=c("Year"), na.rm=TRUE)
D["Age"] <- ifelse(D$variable=="Age2","2", ifelse (D$variable=="Age3","3",
                                                   ifelse (D$variable=="Age4","4",
                                                           ifelse (D$variable=="Age5","5",
                                                                   ifelse (D$variable=="Age6","6",
                                                                           ifelse (D$variable=="Age7","7","8+"))))))

D["EST"]<-D$value
D["Variable"]<-"Estimated"
D<-subset(D, select=c(Year,Age, EST, Variable))

B<-subset(FIGDATA, select=c(Year,obs_sp_comp3, obs_sp_comp4, obs_sp_comp5, obs_sp_comp6, obs_sp_comp7, obs_sp_comp8)) #observed
B["Age3"]<-B$obs_sp_comp3
B["Age4"]<-B$obs_sp_comp4
B["Age5"]<-B$obs_sp_comp5
B["Age6"]<-B$obs_sp_comp6
B["Age7"]<-B$obs_sp_comp7
B["Age8"]<-B$obs_sp_comp8
B["SUM"]<-sum(B$Age3, B$Age4, B$Age5, B$Age6, B$Age7, B$Age8)
B<-subset(B, select=c(Year,Age3, Age4, Age5, Age6, Age7, Age8))
B<- melt(B, id=c("Year"), na.rm=TRUE)
B["Age"] <- ifelse(B$variable=="Age2","2", ifelse (B$variable=="Age3","3",
                                                   ifelse (B$variable=="Age4","4",
                                                           ifelse (B$variable=="Age5","5",
                                                                   ifelse (B$variable=="Age6","6",
                                                                           ifelse (B$variable=="Age7","7","8+"))))))

B["EST"]<-B$value
B["Variable"]<-"Observed"
B<-subset(B, select=c(Year,Age, EST, Variable))
C <- rbind(D,B)
C <- rbind(C,data.frame(Year=2016,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2017,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2018,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2019,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2020,Age=3,EST=NA, Variable="Estimated"))
C <- rbind(C,data.frame(Year=2021,Age=3,EST=NA, Variable="Estimated"))
C$Year<-factor(C$Year,levels=c(
  "1980", "1987","1994","2001", "2008","2015",
  "1981","1988","1995","2002", "2009","2016",
  "1982","1989","1996","2003", "2010", "2017",
  "1983","1990","1997","2004","2011", "2018",
  "1984","1991","1998","2005","2012","2019",
  "1985","1992","1999","2006", "2013","2020",
  "1986","1993","2000","2007","2014","2021"))
g9a<-ggplot(data=C,aes(x=Age, y=EST, fill=Variable))+facet_wrap(~Year,nrow=7, ncol=6)+geom_bar(stat="identity", position="dodge")+
 theme_bw()+xlab ("Age")+ylab("Proportion-at-age")+
  scale_fill_manual(values=c("#E69F00", "#0072B2"))+theme(legend.title=element_blank())+theme(text=element_text(family="Times New Roman", face="bold", size=12))
g9a<-g9a+coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(size=10,colour="black", family="Times New Roman"),
        axis.title.x = element_text(size=14, colour="black", family="Times New Roman"))+
  theme(axis.text.y = element_text(size=10,colour="black", family="Times New Roman"),
        axis.title.y = element_text(size=14,colour="black", family="Times New Roman"))+
  theme(plot.title=element_text(size=rel(1.5),colour="black",vjust =1, family="Times New Roman"))+
  theme(strip.text.x = element_text(size=14,face="bold", family="Times New Roman"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
png(file='Figure 9a.png', res=200, width=9, height=9, units ="in")  
grid.Newpage()
pushViewport(viewport(layout=grid.layout(1,1)))
vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
print(g9a,vp=vplayout(1,1:1)) 
dev.off()
Figure_9a_data<-C
rm(D,B,C,g9a)

#################################################################################################################
#################################################################################################################
#SECTION # III: OUTPUT FILES TO ONE EXCEL WORKBOOK & FORMAT
#################################################################################################################
#################################################################################################################
wb <- createWorkbook()         # create blank workbook
# create different sheets
sheet1 <- createSheet(wb, sheetName="Figure 1") 
sheet2 <- createSheet(wb, sheetName="Figure 2") 
sheet3 <- createSheet(wb, sheetName="Figure 3") 
sheet4 <- createSheet(wb, sheetName="Figure 4")
sheet5 <- createSheet(wb, sheetName="Figure 5")
sheet6 <- createSheet(wb, sheetName="Figure 6")
sheet7 <- createSheet(wb, sheetName="Figure 7") 
sheet8 <- createSheet(wb, sheetName="Figure 8")
sheet9 <- createSheet(wb, sheetName="Figure 9")
sheet10 <- createSheet(wb, sheetName="Figure 10") 
sheet11 <- createSheet(wb, sheetName="Figure 11") 
sheet12 <- createSheet(wb, sheetName="Figure 12") 
sheet13 <- createSheet(wb, sheetName="Figure 13")
sheet14 <- createSheet(wb, sheetName="Figure 14")
sheet15 <- createSheet(wb, sheetName="Figure 15")
sheet16 <- createSheet(wb, sheetName="Figure 8a-Extra")#Sitka only
sheet17 <- createSheet(wb, sheetName="Figure 9a-Extra")#Sitka Only
sheet18 <- createSheet(wb, sheetName="Model Analysis")

setColumnWidth(sheet1, colIndex=17, colWidth=12)
setColumnWidth(sheet1, colIndex=18, colWidth=12)
setColumnWidth(sheet2, colIndex=21, colWidth=15)
setColumnWidth(sheet2, colIndex=22, colWidth=15)
setColumnWidth(sheet3, colIndex=21:26, colWidth=15)
setColumnWidth(sheet4, colIndex=24, colWidth=15)
setColumnWidth(sheet6, colIndex=23:28, colWidth=15)
setColumnWidth(sheet10, colIndex=14, colWidth=15)
setColumnWidth(sheet11, colIndex=12, colWidth=15)
setColumnWidth(sheet12, colIndex=12, colWidth=15)
setColumnWidth(sheet14, colIndex=28, colWidth=25)

rows <- createRow(sheet7,rowIndex=1)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Gear Selectivity")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet7,rowIndex=34)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Maturity")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet7,rowIndex=69)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Survival")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet8,rowIndex=1)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Catch Age Composition")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet9,rowIndex=1)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Spawning Age Composition")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet15,rowIndex=1)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Spawning Age Composition Residuals")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet16,rowIndex=1)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Catch Age Composition")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet17,rowIndex=1)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Spawning Age Composition")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet18,rowIndex=1)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Spawning Age Composition")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

rows <- createRow(sheet18,rowIndex=49)
sheetTitle <- createCell(rows, colIndex=1)
setCellValue(sheetTitle[[1,1]], "Catch Age Composition")
csSheetTitle <- CellStyle(wb) + Font(wb,isBold=TRUE, heightInPoints=12, name="Cambria")
setCellStyle(sheetTitle[[1,1]], csSheetTitle)

cs1 <- CellStyle(wb) + Font(wb, heightInPoints=10, name="Cambria") +
  Alignment(h="ALIGN_CENTER")
addDataFrame(Figure_1_data, sheet1, startRow=1, startColumn=16,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_2_data, sheet2, startRow=1, startColumn=20,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_3_data, sheet3, startRow=1, startColumn=20,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_4_data, sheet4, startRow=1, startColumn=22,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_5_data, sheet5, startRow=1, startColumn=22,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_6_data, sheet6, startRow=1, startColumn=22,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1) 
addDataFrame(Figure_7a_data, sheet7, startRow=1, startColumn=16,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_7b_data, sheet7, startRow=1, startColumn=20,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_7c_data, sheet7, startRow=1, startColumn=24,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_8_data, sheet8, startRow=1, startColumn=20,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_9_data, sheet9, startRow=1, startColumn=16,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1) 
addDataFrame(Figure_10_data, sheet10, startRow=1, startColumn=13,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_11_data, sheet11, startRow=1, startColumn=11,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_12_data, sheet12, startRow=1, startColumn=11,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1) 
addDataFrame(Figure_13_data, sheet13, startRow=1, startColumn=11,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_14_data, sheet14, startRow=1, startColumn=27,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_15_data, sheet15, startRow=1, startColumn=16,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1) 
addDataFrame(Figure_8a_data, sheet16, startRow=1, startColumn=16,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1)
addDataFrame(Figure_9a_data, sheet17, startRow=1, startColumn=16,row.names=FALSE, colnamesStyle=cs1, rownamesStyle=cs1, colStyle=cs1) 

#replace with correct stock
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 1.png", sheet1, startRow=1, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 2.png", sheet2, startRow=1, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 3.png", sheet3, startRow=1, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 4.png", sheet4, startRow=1, startColumn=1)
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 5.png", sheet5, startRow=1, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 6.png", sheet6, startRow=1, startColumn=1)
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 7a.png", sheet7, startRow=2, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 7b.png", sheet7, startRow=36, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 7c.png", sheet7, startRow=71, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 8.png", sheet8, startRow=3, startColumn=1)
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 9.png", sheet9, startRow=3, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 10.png", sheet10, startRow=1, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 11.png", sheet11, startRow=1, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 12.png", sheet12, startRow=1, startColumn=1)
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 13.png", sheet13, startRow=1, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 14.png", sheet14, startRow=1, startColumn=1)
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 15.png", sheet15, startRow=3, startColumn=1) 
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 16.png", sheet18, startRow=2, startColumn=1)
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 17.png", sheet18, startRow=50, startColumn=1) 

#you should get an error here if you are running Craig, Tenakee, or Seymour stocks
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 8a.png", sheet16, startRow=3, startColumn=1)
addPicture("C:/2016 GIT/Sitka/ADMB/211a/Figure 9a.png", sheet17, startRow=3, startColumn=1) 

saveWorkbook(wb, "C:/2016 GIT/Sitka/ADMB/211a/Figures and Data.xlsx")
#remove unneeded files from H drive
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 1.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 2.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 3.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 4.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 5.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 6.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 7c.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 8.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 9.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 10.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 11.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 12.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 13.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 14.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 15.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 7a.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 7b.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 8a.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 9a.png'
if (file.exists(fn)) file.remove(fn)
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 16.png'
if (file.exists(fn)) file.remove(fn)
fn <- 'C:/2016 GIT/Sitka/ADMB/211a/Figure 17.png'
if (file.exists(fn)) file.remove(fn)












































