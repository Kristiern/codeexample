# Examples of plotting 2 best mixtures for stations
# and adding Q-Q-plots and P-P-plots, finding quantiles
# For master's thesis
# Kristi Ernits

library(statmod)
library(actuar)
library(VGAM)
library(truncnorm)
source("functions.R")
source("load_data.R")

# Ristna ----
sub=data2[data2$Parameter=="Wind speed" & data2$Station=="Ristna",]
x=aggregate(Value~Date,FUN=max,data=sub)$Value
datax=as.data.frame(x)
wind(station="Ristna",day=T,estimate=F)
gammaplot(p,shape,rate,max=26,add=T,comps=T)
QQplot2(x,fn=pgammaK,station="Ristna",dist="gamma",max=25)
FF=pgammaK(x,p,shape,rate)
PPplot(x,FF,station="Ristna",dist="gamma")
quantiles(x,pgammaK,max=25)
#
wind(station="Ristna",day=T,estimate=F)
lnormplot(p,meanlog,sdlog,max=26,add=T,comps=T)
QQplot2(x,fn=plnormK,station="Ristna",dist="lognormal",max=25)
FF=plnormK(x,p,meanlog,sdlog)
PPplot(x,FF,station="Ristna",dist="lognormal")
quantiles(x,plnormK,max=25)
