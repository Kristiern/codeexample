# Find summary statistics, make boxplots
# For master's thesis
# Kristi Ernits

source("functions.R")
source("load_data.R")

### Find summmary statistics ----
# Daily maximum of hourly mean (average) wind speeds 
a1=stats(parameter="Wind speed",day=T)
# Hourly mean (average) wind speed
a2=stats(parameter="Wind speed",day=F)
# Daily maximal wind speed
a3=stats(parameter="Max wind speed",day=T)
# Hourly maximal wind speed
a4=stats(parameter="Max wind speed",day=F)
# Daily mean temperature
a5=stats(parameter="Temperature",day=T)
# Hourly mean temperature
a6=stats(parameter="Temperature",day=F)

### Boxplots ----
boxplot2(parameter="Wind speed",day=T)
boxplot2(parameter="Temperature",day=T)
