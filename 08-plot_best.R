# Example of plotting best mixture for station by distribution and overall,
# Q-Q-plots and P-P-plots
# For master's thesis
# Kristi Ernits

library(statmod)
library(actuar)
library(VGAM)
library(truncnorm)
source("functions.R")
source("load_data.R")

# Daily maximum of hourly mean (average) wind speeds ----
wind(station="Ristna",day=T,estimate=F)
legend("topright",col=c(1,2,4,5,6,3,1,2),lty=c(1,2,4:6,6,2,4),lwd=2,cex=0.8,box.lty=0,inset=0.05,
       legend=c("Lognormal mixture","Inverse Gaussian mixture","Gamma mixture",
                "Burr mixture","Inverse Burr mixture","Weibull mixture",
                "Rayleigh mixture","Truncated normal mixture"))
lnormplot(p,meanlog,sdlog,max=26,add=T)
invGaussplot(p,mu,lambda,max=26,add=T)
weibullplot(p,shape,scale,max=26,add=T)
gammaplot(p,shape,rate,max=26,add=T)
burrplot(p,shape1,shape2,scale,max=26,add=T)
invburrplot(p,shape1,shape2,scale,max=26,add=T)
rayleighplot(p,scale,max=26,add=T)
tnorm0plot(p,meann,sdn,max=26,add=T)

# Best mixture with components (wind) ----
wind(station="Ristna",day=T,estimate=F)
gammaplot(p,shape,rate,max=26,add=T,comps=T)
wind(station="Tallinn-Harku",day=T,estimate=F)
tnorm0plot(p,meann,sdn,max=26,add=T,comps=T,short=T)
wind(station="Tõravere",day=T,estimate=F)
gammaplot(p,shape,rate,max=26,add=T,comps=T)
wind(station="Vilsandi",day=T,estimate=F)
invGaussplot(p,mu,lambda,max=26,add=T,comps=T,short=T)
wind(station="Võru",day=T,estimate=F)
gammaplot(p,shape,rate,max=26,add=T,comps=T)

# P-P-plots and Q-Q-plots (wind) ----
# Ristna
sub=data2[data2$Parameter=="Wind speed" & data2$Station=="Ristna",]
x=aggregate(Value~Date,FUN=max,data=sub)$Value
FF=pgammaK(x,p,shape,rate)
PPplot(x,FF,station="Ristna",dist="Gamma")
QQplot2(x,fn=pgammaK,station="Ristna",dist="Gamma",max=25)
 
