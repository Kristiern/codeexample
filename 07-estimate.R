# Find parameter estimates for mixtures
# For master's thesis
# Kristi Ernits

library(flexmix)
library(truncnorm)
source("functions.R")
source("load_data.R")
source("FLXMCdist2.R")

# Daily maximum of hourly mean (average) wind speeds
stations=c("Tallinn-Harku","Ristna","Tõravere","Vilsandi","Võru")
distributions=c("lnorm","invGauss","weibull","gamma","burr","invburr","rayleigh","tnorm0")
for (i in 1:length(stations)){
  for (j in 1:length(distributions)){
    station=stations[i]
    distribution=distributions[j]
    wind(station,day=T,distribution=distribution)
  }
}

# Daily average temperature
temp(station="Ristna",day=T)
temp(station="Tallinn-Harku",day=T)
temp(station="Tõravere",day=T)
temp(station="Vilsandi",day=T)
temp(station="Võru",day=T)

