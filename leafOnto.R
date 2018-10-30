source('readLeafOnto/readLibsLeafOnto.R')
source('readLeafOnto/basicFun.R')
source('readLeafOnto/Vcmax1pt.R')

gammaStar <- 42.75#at 25C according to Bernacchi et al. (2001)
Km <- 710.3203#at 25 C for pO2 210 mmol/mol and Bernacchi et al. (2001)
gasex <- read.csv('dataLeafOnto/spotsINRAgh.csv')
gasex$DateTime <- paste0(gasex$date, ' ', gasex$HHMMSS)
gasex$DateTime <- dmy_hms(as.character(gasex$DateTime))
gasex$Time <- hour(gasex$DateTime)+minute(gasex$DateTime)/60
gasex$rep <- paste0(gasex$round, '_', gasex$spp, '_', gasex$treatment, '_', gasex$leafAge)
gasex$gc <- gasex$Cond/1.6
gasex$CiCa <- gasex$Ci/gasex$CO2S
gasex$vcmax_Tins <- gasex$Photo/(((gasex$Ci-gammaStar)/(gasex$Ci+Km))-0.015)
gasex$vcmax_Tins <- ifelse(gasex$vcmax_Tins<0, NA, gasex$vcmax_Tins)
gasex$vcmax_Tins <- ifelse(gasex$Photo<0, NA, gasex$vcmax_Tins)
gasex$vcmax_Tins <- ifelse(gasex$vcmax_Tins>=100, NA, gasex$vcmax_Tins)
gasex$vcmax_Tins <- ifelse(gasex$spp=='ZEA', NA, gasex$vcmax_Tins)
gasex$vcmax_Tsen <- Vcmax1pt(Photo=gasex$Photo, Ci=gasex$Ci, Tleaf=gasex$Tleaf)
gasex$vcmax_Tsen <- ifelse(gasex$Photo<0, NA, gasex$vcmax_Tsen)
gasex$vcmax_Tsen <- ifelse(gasex$vcmax_Tsen<0, NA, gasex$vcmax_Tsen)
gasex$vcmax_Tsen <- ifelse(gasex$spp=='ZEA', NA, gasex$vcmax_Tsen)
gasex$vcmax_Tsen <- ifelse(gasex$vcmax_Tsen>=500, NA, gasex$vcmax_Tsen)
gasex$Dmmol <- gasex$VpdL*1000/101.6
identify <- unique(gasex[,c('rep','round','spp','treatment','leafAge')])
gasexSumm <- as.data.frame(summarise(group_by(gasex, rep), Dmean=mean.na(VpdL),Dse=s.err(VpdL),
                          decHour=mean.na(Time), DmmolMean=mean.na(Dmmol),DmmolSE=s.err(Dmmol),
                          CondMean=mean.na(Cond), CondSE=s.err(Cond),
                          PhotoMean=mean.na(Photo),PhotoSE=s.err(Photo),
                          CiCaMean=mean.na(CiCa), CiCaSE=s.err(CiCa),
                          gcMean=mean.na(gc), gcSE=s.err(gc),
                          vcmax_TinsMean=mean.na(vcmax_Tins),vcmax_TinsSE=s.err(vcmax_Tins),
                          vcmax_TsenMean=mean.na(vcmax_Tsen),vcmax_TsenSE=s.err(vcmax_Tsen),sampleNo=noNAlength(Photo)))
gasexSumm <- orderBy(~spp, merge(gasexSumm, identify, by='rep', all.x=T, all.y=F))
gasexSumm <- orderBy(~leafAge, gasexSumm)
write.csv(gasexSumm, file='outputLeafOnto/figure1.csv', row.names=F)
gasex$CAPindex <- gasex$Photo/(sqrt(gasex$Dmmol)*(gasex$CO2S-gammaStar))
gasex$MESindex <- sqrt(gasex$Photo/(1.6*gasex$Dmmol*(gasex$CO2S-gammaStar)))
gasex$fround <- as.factor(gasex$roundNum)
jun <- subset(gasex, campaign== "Jun-15")
aug <- subset(gasex, campaign=="Aug-15")
aug <- droplevels(aug, exclude='well-watered')

keep <- read.csv('dataLeafOnto/keepLocalMax.csv')
gasexShort <- subset(merge(gasex, keep[,c('rep','keepMax')], by='rep', all=T), keepMax=='yes')
gasexShort <- subset(gasexShort, treatment=='control' | treatment=='drought')
gasexShortSumm <- as.data.frame(summarise(group_by(gasexShort, rep), vcmax_TsenMean=mean.na(vcmax_Tsen),
                                          vcmax_TsenSE=s.err(vcmax_Tsen)))
gasexShortSumm <- merge(gasexShortSumm, identify, by='rep', all.x=T, all.y=F)

#old code
windows(12,8)
par(mfrow=c(3,2))
with(subset(jun, leafAge=='mature'),boxplot(Photo~spp, main= 'Photo mature Jun'))
with(subset(jun, leafAge=='young'),boxplot(Photo~spp, main= 'Photo young Jun'))
with(subset(aug, leafAge=='mature' & treatment=='control'),
     boxplot(Photo~spp, main= 'Photo mature Aug Cont'))
with(subset(aug, leafAge=='young' & treatment=='control'),
     boxplot(Photo~spp, main= 'Photo young Aug Cont'))
with(subset(aug, leafAge=='mature' & treatment=='drought'),
     boxplot(Photo~spp, main= 'Photo mature Aug Cont'))
with(subset(aug, leafAge=='young' & treatment=='drought'),
     boxplot(Photo~spp, main= 'Photo young Aug Cont'))

windows(12,8)
par(mfrow=c(3,2))
with(subset(jun, leafAge=='mature'),boxplot(Cond~spp, main= 'Cond mature Jun'))
with(subset(jun, leafAge=='young'),boxplot(Cond~spp, main= 'Cond young Jun'))
with(subset(aug, leafAge=='mature' & treatment=='control'),
     boxplot(Cond~spp, main= 'Cond mature Aug Control'))
with(subset(aug, leafAge=='young' & treatment=='control'),
     boxplot(Cond~spp, main= 'Cond young Aug Control'))
with(subset(aug, leafAge=='mature' & treatment=='drought'),
     boxplot(Cond~spp, main= 'Cond mature Aug Drought'))
with(subset(aug, leafAge=='young' & treatment=='drought'),
     boxplot(Cond~spp, main= 'Cond young Aug Drought'))

windows(12,8)
par(mfrow=c(2,2))
hist(jun$Cond)
hist(jun$Photo)
hist(aug$Cond)
hist(aug$Photo)

anova(lm(Photo~spp*leafAge, data=jun))
windows(12,8)
par(mfrow=c(2,2))
plot(lm(Photo~spp*leafAge, data=jun))
anova(lm(Photo~spp*treatment*leafAge, data=aug))
windows(12,8)
par(mfrow=c(2,2))
plot(lm(Photo~spp*leafAge*treatment, data=aug))

anova(lm(Cond~spp*leafAge, data=jun))
windows(12,8)
par(mfrow=c(2,2))
plot(lm(Cond~spp*leafAge, data=jun))
anova(lm(Cond~spp*treatment*leafAge, data=aug))
windows(12,8)
par(mfrow=c(2,2))
plot(lm(Cond~spp*leafAge*treatment, data=aug))
