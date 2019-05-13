# Read in data
# For master's thesis
# Kristi Ernits

library(openxlsx)
library(tidyr)

# Read in data
raw1=read.xlsx("file.xlsx",
              sheet=1,startRow=5,colNames=FALSE,detectDates=FALSE)
raw2=read.xlsx("file.xlsx",
              sheet=2,startRow=3,colNames=FALSE,detectDates=FALSE)
raw=rbind(raw1[,-23],raw2)
colnames(raw)=c("Year","Month","Day","Time",
                "TH_TEMP","TH_WIND","TH_MWIN",
                "NJ_TEMP","NJ_WIND","NJ_MWIN",
                "RI_TEMP","RI_WIND","RI_MWIN",
                "TO_TEMP","TO_WIND","TO_MWIN",
                "VI_TEMP","VI_WIND","VI_MWIN",
                "VO_TEMP","VO_WIND","VO_MWIN")
# Remove duplicates
raw=unique(raw)
# Correct time
raw$Hour=24*raw$Time 
# Add date
raw$Date=paste(raw$Year,raw$Month,raw$Day,sep="-")
# Transpose
raw_t=gather(raw,key=Varname,value=Value,
             TH_TEMP,TH_WIND,TH_MWIN,
             NJ_TEMP,NJ_WIND,NJ_MWIN,
             RI_TEMP,RI_WIND,RI_MWIN,
             TO_TEMP,TO_WIND,TO_MWIN,
             VI_TEMP,VI_WIND,VI_MWIN,
             VO_TEMP,VO_WIND,VO_MWIN)
# Add station and parameter
raw_t$Station=factor(substr(raw_t$Varname,1,2),
                     levels=c("TH","NJ","RI","TO","VI","VO","N"),
                     labels=c("Tallinn-Harku","Narva-Jõesuu","Ristna",
                              "Tõravere","Vilsandi","Võru","Narva"))
raw_t$Parameter=factor(substr(raw_t$Varname,4,7),
                       levels=c("TEMP","WIND","MWIN"),
                       labels=c("Temperature","Wind speed","Max wind speed"))
# Narva-Jõesuu from 19.12.2013 12:00 UTC in Narva
raw_t$Station[raw_t$Station=="Narva-Jõesuu" & 
              ((raw_t$Year>2013) | 
               (raw_t$Year==2013 & raw_t$Month==12 & raw_t$Day>19) |
               (raw_t$Year==2013 & raw_t$Month==12 & raw_t$Day==19 & raw_t$Hour>11))]="Narva"
# Final data
data=raw_t[,c("Year","Date","Month","Day","Hour","Station","Parameter","Value")]
save(data,file="data.RData")

