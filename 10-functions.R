# Functions
# For master's thesis
# Kristi Ernits

###
### Mixiture characteristics -----
###
mix.mean <- function(p,means){
  mixmean=sum(p*means)
  return(mixmean)
}

mix.var <- function(p,means,vars){
  mixvar=sum(p*(means**2+vars))-(mix.mean(p,means))**2
  return(mixvar)
}

mix.skew <- function(p,means,vars,m3s){
  mixmean=mix.mean(p,means)
  mixvar=mix.var(p,means,vars)
  mixskew=sum(pi*((means-mixmean)**3+3*(means-mixmean)*vars+m3s))/(sqrt(mixvar))**3
  return(mixskew)
}

###
### Cumulative distribution functions for mixtures ----
###
pburrK <- function(q,px=p,shape1x=shape1,shape2x=shape2,scalex=scale,a=0){
  K=length(px)
  FF=q*0
  for (k in 1:K){
    Fk=px[k]*pburr(q,shape1=shape1x[k],shape2=shape2x[k],scale=scalex[k])
    FF=FF+Fk}
  return(FF-a)
}

pgammaK <- function(q,px=p,shapex=shape,ratex=rate,a=0){
  K=length(px)
  FF=q*0
  for (k in 1:K){
    Fk=px[k]*pgamma(q,shape=shapex[k],rate=ratex[k])
    FF=FF+Fk}
  return(FF-a)
}

pinvGaussK <- function(q,px=p,mux=mu,lambdax=lambda,a=0){
  K=length(px)
  FF=q*0
  for (k in 1:K){
    Fk=px[k]*pinvgauss(q,mean=mux[k],shape=lambdax[k])
    FF=FF+Fk}
  return(FF-a)
}

pnormK <- function(q,px=p,meannx=meann,sdnx=sdn,a=0){
  K=length(px)
  FF=q*0
  for (k in 1:K){
    Fk=px[k]*pnorm(q,mean=meannx[k],sd=sdnx[k])
    FF=FF+Fk}
  return(FF-a)
}

ptnorm0K <- function(q,px=p,meannx=meann,sdnx=sdn,a=0){
  K=length(px)
  FF=q*0
  for (k in 1:K){
    Fk=px[k]*ptruncnorm(q,a=0,b=Inf,mean=meannx[k],sd=sdnx[k])
    FF=FF+Fk}
  return(FF-a)
}

plnormK <- function(q,px=p,meanlogx=meanlog,sdlogx=sdlog,a=0){
  K=length(px)
  FF=q*0
  for (k in 1:K){
    Fk=px[k]*plnorm(q,meanlog=meanlogx[k],sdlog=sdlogx[k])
    FF=FF+Fk}
  return(FF-a)
}

prayleighK <- function(q,px=p,scalex=scale,a=0){
  K=length(px)
  FF=q*0
  for (k in 1:K){
    Fk=px[k]*prayleigh(q,scale=scalex[k])
    FF=FF+Fk}
  return(FF-a)
}

###
### Find quantiles ----
###
quantiles <- function(x,fn,min=0,max=15){
  theoretical=c(uniroot(fn,interval=c(min,max),a=0.01)$root,
                uniroot(fn,interval=c(min,max),a=0.05)$root,
                uniroot(fn,interval=c(min,max),a=0.25)$root,
                uniroot(fn,interval=c(min,max),a=0.50)$root,
                uniroot(fn,interval=c(min,max),a=0.75)$root,
                uniroot(fn,interval=c(min,max),a=0.95)$root, 
                uniroot(fn,interval=c(min,max),a=0.99)$root)
  empirical=quantile(x,probs=c(0.01,0.05,0.25,0.50,0.75,0.95,0.99))
  relativedif=100*abs(theoretical-empirical)/abs(empirical)
  tog=rbind(theoretical,empirical,relativedif)
  print(tog)
}

###
### Functions for histograms ----
###
windhist <- function(x,station,tit){ # Mean wind speed
  h=hist(x,breaks=seq(0,26,by=0.5),plot=F)
  ymax=(ceiling(100*max(h$density))+ceiling(100*max(h$density))%%2)/100+0.02
  hist(x,breaks=seq(0,26,by=0.5),freq=F,
       ylim=c(0,ymax),axes=F,border="grey",
       xlab="",ylab="Wind speed frequency",main="")
  title(xlab="Wind speed (m/s)",line=2)
  axis(1,at=seq(0,26,by=4))
  axis(2,at=seq(0,ymax,by=0.02),las=1)
  title(paste("Station:",station),font.main=1,cex.main=1,line=0,adj=0.9)
  title(paste(tit,"mean wind speed"),font.main=1,cex.main=1,line=1,adj=0.9)
}

maxwindhist <- function(x,station,tit){ # Max wind speed
  h=hist(x,breaks=seq(0,34,by=0.5),plot=F)
  ymax=(ceiling(100*max(h$density))+ceiling(100*max(h$density))%%2)/100+0.02
  hist(x,breaks=seq(0,34,by=0.5),freq=F,
       ylim=c(0,ymax),axes=F,border="grey",
       xlab="",ylab="Wind speed frequency",main="")
  title(xlab="Wind speed (m/s)",line=2)
  axis(1,at=seq(0,34,by=4))
  axis(2,at=seq(0,ymax,by=0.02),las=1)
  title(paste("Station:",station),font.main=1,cex.main=1,line=0,adj=0.9)
  title(paste(tit,"maximal wind speed"),font.main=1,cex.main=1,line=1,adj=0.9)
}

temphist <- function(x,station,tit){ # Mean temperature
  h=hist(x,breaks=seq(-35,35,by=1),plot=F)
  ymax=(ceiling(100*max(h$density))+ceiling(100*max(h$density))%%2)/100+0.02
  hist(x,breaks=seq(-35,35,by=1),freq=FALSE,
       ylim=c(0,ymax),axes=F,border="grey",
       xlab="",ylab="Temperature frequency",main="")
  title(xlab="Temperature (°C)",line=2)
  axis(1,at=seq(-35,35,by=10))
  axis(2,at=seq(0,ymax,by=0.02),las=1)
  title(paste("Station:",station),font.main=1,cex.main=1,line=0,adj=0.1)
  title(paste(tit,"mean temperature"),font.main=1,cex.main=1,line=1,adj=0.1)
}

###
### Boxplots ----
###
boxplot2 <- function(parameter="Wind speed",day=F){
  sub=data2[data2$Parameter==parameter,]
  if (day){
    if (parameter=="Temperature"){ 
      x=aggregate(Value~Date+Station,FUN=mean,data=sub)
      tit="Daily"}
    else { 
      x=aggregate(Value~Date+Station,FUN=max,data=sub)
      tit="Daily maximal hourly"}}
  else {x=sub; tit="Hourly"}
  x=x[!(x$Station %in% c("Narva","Narva-Jõesuu")),]
  if (parameter=="Wind speed"){ unit="(m/s)"}
  else {unit="(°C)"}
  boxplot(Value~as.character(Station),data=x,horizontal=F,
          ylab=paste(parameter,unit),main="",xlab="")
  title(paste(tit,"mean",tolower(parameter)),font.main=1,cex.main=1,line=1,adj=0.9)
  points(aggregate(x$Value,by=list(x$Station),FUN=mean)$x,pch=5,col=1)
}

###
### Function for fitting  ----
###
est <- function(datax,max.k,nrep,distribution){
  set.seed(123)
  print("--------------------------")
  print(distribution)
  if (distribution=="rayleigh"){
    fit=stepFlexmix(x~1,data=datax,k=1:max.k,model=FLXMCdist2(dist=distribution),
                    control=list(iter=1000,minprior=0.01),nrep=nrep)}
  else if (distribution=="tnorm0"){
    fit=stepFlexmix(x~1,data=datax,k=1:max.k,model=FLXMCdist2(dist=distribution,method="BFGS"),
                    control=list(iter=1000,minprior=0.01),nrep=nrep)}
  else {
    fit=stepFlexmix(x~1,data=datax,k=1:max.k,model=FLXMCdist1(dist=distribution,method="BFGS"),
                    control=list(iter=1000,minprior=0.01),nrep=nrep)}
  print(fit)
  best=getModel(fit,which="BIC")
  print("Prior")
  print(prior(best))
  print("Parameter estimates")
  print(parameters(best))
}

###
### Function to find summary statistics ----
###
stats <- function(parameter="Wind speed",day=F){
  sub=data2[data2$Parameter==parameter,]
  if (day){
    if (parameter=="Temperature"){ x=aggregate(Value~Date+Station,FUN=mean,data=sub)}
    else { x=aggregate(Value~Date+Station,FUN=max,data=sub)}}
  else {x=sub}
  a=tapply(x$Value,x$Station,summary)
  tab=c()
  for (i in 1:length(a)){
    tab=rbind(tab,a[[i]])}
  n=tapply(x$Value,x$Station,function(x) length(x))
  if (day){ nmiss=5479-n }
  else { nmiss=131496-n }
  tab=cbind(n,nmiss,tab,tapply(x$Value,x$Station,sd))
  rownames(tab)=names(a)
  colnames(tab)=c("n","nmiss","Minimum","Q1","Median","Mean","Q3","Maximum","SD")
  return(tab)
}

###
### Functions for extracting data + histogram + estimate  ----
###
wind <- function(station="Ristna",day=F,estimate=T,max.k=6,nrep=10,distribution="lnorm"){
  # Data
  sub=data2[data2$Parameter=="Wind speed" & data2$Station==station,]
  if (day){ x=aggregate(Value~Date,FUN=max,data=sub)$Value; tit="Daily maximal hourly"}
  else { x=sub$Value; tit="Hourly"}
  x=replace(x,x==0,0.05) 
  datax=as.data.frame(x)
  # Histogram
  windhist(x,station,tit)
  # Estimate
  if (estimate) {est(datax,max.k=max.k,nrep=nrep,distribution=distribution)}
}

maxwind <- function(station="Ristna",day=F,estimate=T,max.k=8,nrep=10,distribution="lnorm"){ 
  # Data
  sub=data2[data2$Parameter=="Max wind speed" & data2$Station==station,]
  if (day){ x=aggregate(Value~Date,FUN=max,data=sub)$Value; tit="Daily"}
  else { x=sub$Value; tit="Hourly"}
  x=replace(x,x==0,0.05) 
  datax=as.data.frame(x)
  # Histogram
  maxwindhist(x,station,tit)
  # Estimate
  #if (estimate) {est(datax,max.k=max.k,nrep=nrep,distribution=distribution)}
}

temp <- function(station="Ristna",day=F,estimate=T,max.k=8,nrep=10){
  # Data
  sub=data2[data2$Parameter=="Temperature" & data2$Station==station,]
  if (day){ x=aggregate(Value~Date,FUN=mean,data=sub)$Value; tit="Daily"}
  else { x=sub$Value; tit="Hourly"}
  datax=as.data.frame(x)
  # Histogram
  temphist(x,station,tit)
  # Estimate
  if (estimate) {
      set.seed(123)
      print("--------------------------")
      print("Normal")
      fit=stepFlexmix(x~1,data=datax,k=1:max.k,
                      model=FLXMCnorm1(),
                      control=list(iter=1000,minprior=0.01),nrep=nrep)
      print(fit)
      best=getModel(fit,which="BIC")
      print("Prior")
      print(prior(best))
      print("Parameter estimates")
      print(parameters(best))}
}

###
### Functions to plot mixtures or add mixtures to histogram ----
###
lnormplot <- function(p,meanlog,sdlog,add=T,max=26,max2=1,comps=F,short=F){
  K=length(p)
  xx=seq(0,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(0,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*dlnorm(xx,meanlog=meanlog[k],sd=sdlog[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=1,lty=1,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    if (short){
    legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
           legend=c(expression(paste("Density of lognormal mixture, f(y|",Phi,")")),
                    paste("Mix. comp., ",round(p,2),
                          "g(y|",round(meanlog,2),",",round(sdlog,2),")",sep="")))}
    else{
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
             legend=c(expression(paste("Density of lognormal mixture, f(y|",Phi,")")),
                      paste("Mixture component, ",round(p,2),
                            "g(y|",round(meanlog,2),",",round(sdlog,2),")",sep="")))}}
  means=exp(meanlog+sdlog**2/2)
  vars=exp(2*meanlog+sdlog**2)*(exp(sdlog**2)-1)
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)
}

invGaussplot <- function(p,mu,lambda,add=T,max=26,max2=1,comps=F,short=F){
  K=length(p)
  xx=seq(0,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(0,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*dinvgauss(xx,mean=mu[k],shape=lambda[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=2,lty=2,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    if (short){
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
           legend=c(expression(paste("Density of inv. Gaussian mix., f(y|",Phi,")")),
                    paste("Mix. comp., ",round(p,2),
                          "g(y|",round(mu,2),",",round(lambda,2),")",sep="")))}
    else {
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
             legend=c(expression(paste("Density of inverse Gaussian mixture, f(y|",Phi,")")),
                      paste("Mixture component, ",round(p,2),
                            "g(y|",round(mu,2),",",round(lambda,2),")",sep="")))}}
  means=mu
  vars=mu**3/lambda
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)
  m3s=3*mu**5/lambda**2
  print(paste("Skewness: ",mix.skew(p,means,vars,m3s))); print(m3s)}

weibullplot <- function(p,shape,scale,add=T,max=26,max2=1,comps=F){
  K=length(p)
  xx=seq(0,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(0,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*dweibull(xx,shape=shape[k],scale=scale[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=3,lty=6,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
           legend=c(expression(paste("Density of Weibull mixture, f(y|",Phi,")")),
                    paste("Mixture component, ",round(p,2),
                          "g(y|",round(shape,2),",",round(scale,2),")",sep="")))}
  means=scale*gamma(1+1/shape)
  vars=scale**2 *gamma(1+2/shape) - (scale*gamma(1+1/shape))**2
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)
}

gammaplot <- function(p,shape,rate,add=T,max=26,max2=1,comps=F){
  K=length(p)
  xx=seq(0,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(0,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*dgamma(xx,shape=shape[k],rate=rate[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=4,lty=4,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
           legend=c(expression(paste("Density of gamma mixture, f(y|",Phi,")")),
                    paste("Mixture component, ",round(p,2),
                          "g(y|",round(shape,2),",",round(rate,2),")",sep="")))}
  means=shape/rate
  vars=shape/rate**2
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)
}

burrplot <- function(p,shape1,shape2,scale,add=T,max=26,max2=1,comps=F,short=F){
  K=length(p)
  xx=seq(0,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(0,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*dburr(xx,shape1=shape1[k],shape2=shape2[k],scale=scale[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=5,lty=5,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    if (short){
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
             legend=c(expression(paste("Density of Burr mixture, f(y|",Phi,")")),
                      paste("Mix. comp., ",round(p,2),
                            "g(y|",round(shape1,2),",",round(shape2,2),",",round(scale,2),")",
                            sep="")))}
    else {
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
             legend=c(expression(paste("Density of Burr mixture, f(y|",Phi,")")),
                      paste("Mixture component, ",round(p,2),
                            "g(y|",round(shape1,2),",",round(shape2,2),",",round(scale,2),")",
                            sep="")))}}
  if (all(shape1*shape2>1)){
  means=scale*shape1*beta(shape1-1/shape2,1+1/shape2)
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  if (all(shape1*shape2>2)){
  vars=scale**2 * shape1 *(beta(shape1-2/shape2,1+2/shape2)-shape1*(beta(shape1-1/shape2,1+1/shape2))**2)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)}}
}

invburrplot <- function(p,shape1,shape2,scale,add=T,max=26,max2=1,comps=F,short=F){
  K=length(p)
  xx=seq(0,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(0,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*dinvburr(xx,shape1=shape1[k],shape2=shape2[k],scale=scale[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=6,lty=6,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    if (short){
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
             legend=c(expression(paste("Density of inverse Burr mixture, f(y|",Phi,")")),
                      paste("Mix. comp., ",round(p,2),
                            "g(y|",round(shape1,2),",",round(shape2,2),",",round(scale,2),")",
                            sep="")))}
    else {
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
             legend=c(expression(paste("Density of inverse Burr mixture, f(y|",Phi,")")),
                      paste("Mixture component., ",round(p,2),
                            "g(y|",round(shape1,2),",",round(shape2,2),",",round(scale,2),")",
                            sep="")))}}
  if (all(shape2>1)){
  means=scale*shape1*beta(shape1+1/shape2,1-1/shape2)
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  if (all(shape2>2)){
  vars=scale**2 * shape1 * (beta(shape1+2/shape2,1-2/shape2)-shape1*(beta(shape1+1/shape2,1-1/shape2))**2)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)}}
}

rayleighplot <- function(p,scale,add=T,max=26,max2=1,comps=F){
  K=length(p)
  xx=seq(0,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(0,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*drayleigh(xx,scale=scale[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=1,lty=2,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
           legend=c(expression(paste("Density of Rayleigh mixture, f(y|",Phi,")")),
                    paste("Mixture component, ",round(p,2),
                          "g(y|",round(scale,2),")",sep="")))}
  means=sqrt(pi/2)*scale
  vars=(2-pi/2)*scale**2
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)
}

normplot <- function(p,meann,sdn,add=T,comps=F,max=35,max2=1,min=-35){
  K=length(p)
  xx=seq(min,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(min,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*dnorm(xx,mean=meann[k],sd=sdn[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=1,lty=1,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
           legend=c(expression(paste("Density of normal mixture, f(y|",Phi,")")),
                    paste("Mixture component, ",round(p,2),
                          "g(y|",round(meann,2),",",round(sdn,2),")",sep="")))}
  means=meann
  vars=sdn**2
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)
}

tnorm0plot <- function(p,meann,sdn,add=T,comps=F,max=26,max2=1,short=F){
  K=length(p)
  xx=seq(0,max,by=0.001)
  f=0*xx
  if (!add) {
    plot(0,type='n',xlim=c(0,max),ylim=c(0,max2),xlab="y",xaxs="i",
         ylab=expression(paste("Mixture density, f(y|",Phi,")")),cex.lab=0.8,cex.axis=0.8)}
  for (k in 1:K){
    fk=p[k]*dtruncnorm(xx,a=0,b=Inf,mean=meann[k],sd=sdn[k])
    f=f+fk
    if (!add) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}
    if (comps) {lines(xx,fk,type='l',lty=2,lwd=1,col=k+1)}}
  if (add & !comps) {lines(xx,f,col=2,lty=4,lwd=1)}
  else {
    lines(xx,f,col=1,lty=1,lwd=1)
    if (short){
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
             legend=c(expression(paste("Density of trunc. normal mix., f(y|",Phi,")")),
                      paste("Mix. comp., ",round(p,2),
                            "g(y|",round(meann,2),",",round(sdn,2),")",sep="")))}
    else {
      legend("topright",col=1:(K+1),lty=c(1,rep(2,K)),cex=0.8,inset=0.05,bty='n',xpd=NA,
             legend=c(expression(paste("Density of truncated normal mixture, f(y|",Phi,")")),
                      paste("Mixture component, ",round(p,2),
                            "g(y|",round(meann,2),",",round(sdn,2),")",sep="")))}}
  means=meann+sdn/sqrt(2*pi)/pnorm(meann/sdn)*exp(-meann**2/2/sdn**2)
  vars=sdn**2-meann*sdn/sqrt(2*pi)/pnorm(meann/sdn)*exp(-meann**2/2/sdn**2)-sdn**2/(2*pi)/((pnorm(meann/sdn))**2)*exp(-meann**2/sdn**2)
  print(paste("Mean: ", mix.mean(p,means))); print(means)
  print(paste("Variance: ",mix.var(p,means,vars))); print(vars)
}

###
### P-P-plot ----
###
PPplot <- function(x,FF,station,dist,temp=F,title=T){
  Fe=ecdf(x)
  plot(Fe(x),FF,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",
       panel.first=c(abline(0,1,lty=2,col=2)))
  title(xlab="Observed",ylab=paste("Expected (",dist," mixture)",sep=""),line=2)
  if (title){
  title(paste("Station:",station),font.main=1,cex.main=1,line=1,adj=0.9)
  if (temp) {title("Daily mean temperature",font.main=1,cex.main=1,line=2,adj=0.9)}
  else {title("Daily maximal hourly mean wind speed",font.main=1,cex.main=1,line=2,adj=0.9)}}
}

###
### Q-Q-plot ----
###
QQplot2 <- function(x,fn,station,dist,min=0,max=15,temp=F,title=T){
  n=length(x)
  tqs=rep(0,n)
  for (i in 1:n){
    a=(i-0.5)/n
    tqs[i]=uniroot(fn,interval=c(min,max),a=a)$root}
  plot(sort(x),tqs,xlim=c(min,max),ylim=c(min,max),xlab="",ylab="",
       panel.first=c(abline(0,1,lty=2,col=2)))
  title(xlab="Empirical quantiles",ylab=paste("Theoretical quantiles (",dist," mixture)",sep=""),line=2)
  if (title) {
  title(paste("Station:",station),font.main=1,cex.main=1,line=1,adj=0.9)
  if (temp) {title("Daily mean temperature",font.main=1,cex.main=1,line=2,adj=0.9)}
  else {title("Daily maximal hourly mean wind speed",font.main=1,cex.main=1,line=2,adj=0.9)}}
}
