meteoGH_pre <- read.csv('dataLeafOnto/climate_glasshouse_pre.csv')
dates <- read.csv('dataLeafOnto/datesDOY.csv')
meteoGH_pre$Time <- meteoGH_pre$hhmm
meteoGH_pre[which(meteoGH_pre$hhmm < 1000), 'Time'] <- c(paste0("0", meteoGH_pre[which(meteoGH_pre$hhmm < 1000), 'hhmm']))
meteoGH_pre[which(meteoGH_pre$hhmm < 100), 'Time'] <- c(paste0("00", meteoGH_pre[which(meteoGH_pre$hhmm < 100), 'hhmm']))
meteoGH_pre[which(meteoGH_pre$hhmm == 2400), 'Time'] <- "0000"
meteoGH_pre[which(meteoGH_pre$hhmm == 2400), 'DOY'] <- meteoGH_pre[which(meteoGH_pre$hhmm == 2400), 'DOY']+1
meteoGH_pre <- merge(meteoGH_pre, dates, by = 'DOY', all.x=TRUE, all.y=FALSE)
meteoGH_pre$DateTime <- paste0(meteoGH_pre[, 'Date'], " ", meteoGH_pre[, 'Time'])
meteoGH_pre$DateTime <- dmy_hm(as.character(meteoGH_pre$DateTime))
#data are OK only from Thurs 4 June 2015, DOY = 155
meteoGH_pre <- subset(meteoGH_pre, DOY >= 155)
#PAR sensor in place only from Fri 5 Jun 2015, DOY = 156
meteoGH_pre[which(meteoGH_pre$DOY <= 155), "Par_Den_AVG"] <- NA
meteoGH_pre$RH <- c(rep(NA, times=nrow(meteoGH_pre)))
meteoGH_pre <- meteoGH_pre[, c('DateTime', 'DOY', "T_C_Avg", "RH", "Par_Den_AVG")]
names(meteoGH_pre) <- c("DateTime", "DOY", "Temp", "RH", "PAR")

meteoGH <- read.csv('dataLeafOnto/climate_glasshouse.csv')
meteoGH$Time <- meteoGH$hhmm
meteoGH[which(meteoGH$hhmm < 1000), 'Time'] <- c(paste0("0", meteoGH[which(meteoGH$hhmm < 1000), 'hhmm']))
meteoGH[which(meteoGH$hhmm < 100), 'Time'] <- c(paste0("00", meteoGH[which(meteoGH$hhmm < 100), 'hhmm']))
meteoGH <- merge(meteoGH, dates, by='DOY',all.x=TRUE, all.y=FALSE)
meteoGH$DateTime <- paste0(meteoGH[, 'Date'], " ", meteoGH[,  'Time'])
meteoGH$DateTime <- dmy_hm(as.character(meteoGH$DateTime))
meteoGH <- meteoGH[, c('DateTime', "DOY", 'T_C_Avg', 'RH_frac_Avg', 'PAR_Den_Avg')]
names(meteoGH) <- c("DateTime", "DOY", "Temp", "RH", "PAR")

meteoGH <- rbind(meteoGH_pre, meteoGH)
meteoGH$Date <- as.Date(meteoGH$DateTime)
meteoGH$VPD <- calcVPD(meteoGH$RH, meteoGH$Temp)
meteoCampaignJun <- subset(meteoGH, Date==as.Date('2015-06-18') | Date==as.Date('2015-06-22') | 
                             Date==as.Date('2015-06-24'))
meteoCampaignAug <- subset(meteoGH, Date==as.Date('2015-08-05') | Date==as.Date('2015-08-06') | 
                             Date==as.Date('2015-08-20') | Date==as.Date('2015-08-24'))

meteoExt <- read.csv('dataLeafOnto/climate_ext.csv')
meteoExt$DateTime <- paste0(meteoExt$Date, " ", meteoExt$Hour)
meteoExt$DateTime <- nearestTimeStep(dmy_hm(as.character(meteoExt$DateTime)),10)
meteo <- merge(meteoGH, meteoExt, by='DateTime',all=TRUE)

meteoRad <- read.csv('dataLeafOnto/climate_radio.csv')
meteoRad$DateTime <- dmy_hm(meteoRad$DateTime)
meteo <- merge(meteo, meteoRad[, c('DateTime','Temp_radio','RH_radio','PAR_radio')],by='DateTime',all=TRUE)

meteoElsa <- read.csv('dataLeafOnto/climate_int_Elsa.csv')
meteoElsa$DateTime <- nearestTimeStep(dmy_hm(as.character(paste0(meteoElsa$Date, " ", meteoElsa$Hour))),10)
meteo <- merge(meteo, meteoElsa[,c('DateTime','Temp_elsa',"RH_elsa")], by='DateTime', all=TRUE)
meteo$Date <- as.Date(meteo$DateTime)
meteo$PARconv <- meteo$PAR * 10 * 60 /1000000
meteo$radConv <- meteo$SolarRad_ext *10 /1000000
tempModel <- lm(Temp~Temp_ext, data=meteo)
meteo$Tpred <- coefficients(tempModel)[1]+coefficients(tempModel)[2]*meteo$Temp_ext
meteoDay <- subset(meteo, SolarRad_ext>0)
meteoNight <- subset(meteo, SolarRad_ext==0)

meteoDOY <- as.data.frame(summarize(group_by(meteo,Date),TempMean=mean(Temp),TempMax=max(Temp),Temp=min(Temp),VPD=mean(VPD),
                                    PAR=sum(PARconv),TextMean=mean(Temp_ext),TextMax=max(Temp_ext),TextMin=min(Temp_ext),
                                    radExt=sum(radConv),TpredMean=mean(Tpred),
                                    TpredMax=max(Tpred),TpredMin=min(Tpred)))
dayPAR <- subset(meteo, PAR>=20)
summPAR <- summaryBy(PAR~Date, FUN=mean.na, data=meteo)
PARgh <-c(mean(summPAR$PAR.mean.na), s.err(summPAR$PAR.mean.na), length(summPAR$PAR.mean.na))

windows(12,10)
par(mfrow=c(2,2), mar=c(2,4,2,2))
with(meteoGH, plot(Temp~DateTime, type="l", col= "red", xlab=" ", ylab="T (C)"))
with(meteoGH, plot(RH~DateTime, type="l", col= "blue", xlab=" ", ylab="RH (%)"))
with(meteoGH, plot(VPD~DateTime, type="l", col= "forestgreen", xlab=" ", ylab="D (kPa)"))
with(meteoGH, plot(PAR~DateTime, type="l", col= "black", xlab=" ", ylab="PAR"))