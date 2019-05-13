# Examples
# For master's thesis
# Kristi Ernits

library(statmod)
library(actuar)
source("functions.R")

# Lognormal
p=c(0.4,0.3,0.3)
meanlog=c(1.5,2.5,3)
sdlog=c(0.25,0.25,0.5)
lnormplot(p,meanlog,sdlog,add=F,max=40,max2=0.15)
# Inverse Gaussian
p=c(0.6,0.4)
mu=c(6,1)
lambda=c(2,10)
invGaussplot(p,mu,lambda,add=F,max=10,max2=0.8)
# Gamma
p=rep(0.25,4)
shape=c(2,4,8,16) #alpha
rate=rep(2,4) #lambda
gammaplot(p,shape,rate,add=F,max=15,max2=0.25)
# Burr
p=c(0.8,0.2)
shape1=c(4,2) # alpha
shape2=c(1,4) # gamma
scale=c(2,3) # theta
burrplot(p,shape1,shape2,scale,add=F,max=6,max2=1.6)
# Inverse Burr
p=c(0.3,0.3,0.4)
shape1=c(0.5,2,2) # alpha
shape2=c(1,3,4) # gamma
scale=c(5,1,4) # theta
invburrplot(p,shape1,shape2,scale,add=F,max=10,max2=0.8)
# Weibull
p=c(0.8,0.2)
scale=c(3,1)
shape=c(2,3)
weibullplot(p,shape,scale,add=F,max=8,max2=0.4)
# Normal
set.seed(12)
x=rnorm(200,0,15)
hist(x,breaks=seq(-40,40,by=2.5),freq=FALSE,
     ylim=c(0,0.04),axes=F,border="grey",
     xlab="",ylab="Frequency",main="")
title(xlab="x",line=2)
axis(1,at=seq(-40,40,by=20),cex.axis=1.5)
axis(2,at=seq(0,0.04,by=0.02),las=1)

hist(x,breaks=seq(-40,40,by=5),freq=FALSE,
     ylim=c(0,0.04),axes=F,border="grey",
     xlab="",ylab="Frequency",main="")
title(xlab="x",line=2)
axis(1,at=seq(-40,40,by=20))
axis(2,at=seq(0,0.04,by=0.02),las=1)
