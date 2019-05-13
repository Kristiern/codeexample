# Examples of plotting histograms for stations
# For master's thesis
# Kristi Ernits

source("functions.R")
source("load_data.R")

# Daily maximum of hourly mean (average) wind speeds 
wind(station="Ristna",day=T,estimate=F)
# Hourly mean (average) wind speeds 
wind(station="Ristna",day=F,estimate=F)
# Daily maximum wind speeds
maxwind(station="Ristna",day=T,estimate=F)
# Hourly maximum wind speeds
maxwind(station="Ristna",day=F,estimate=F)
# Daily mean temperature
temp(station="Ristna",day=T,estimate=F)
# Hourly mean temperature
temp(station="Ristna",day=F,estimate=F)
