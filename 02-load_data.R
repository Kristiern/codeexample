# Load data
# For master's thesis
# Kristi Ernits

### Load weather data ----
load("data.RData") # 2366928
missing=data[is.na(data$Value),]
data2=data[!is.na(data$Value),] 

