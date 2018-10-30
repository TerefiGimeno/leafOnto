windows(12,8)
par(mfrow=c(2,2))
boxplot(jun$vcmax_Tsen~(jun$leafAge*jun$spp), ylim=c(0,100))
boxplot(aug$vcmax_Tsen~(aug$treatment*aug$spp), ylim=c(0,100))
boxplot(log(jun$vcmax_Tsen+1)~(jun$leafAge*jun$spp), ylim=c(0,5))
boxplot(log(aug$vcmax_Tsen+1)~(aug$treatment*aug$spp), ylim=c(0,5))

vcmax_Tsen <- as.data.frame(summarise(group_by(gasex, spp, leafAge, treatment), 
                                     vcmax_TsenMean=mean.na(vcmax_Tsen),vcmax_TsenSE=s.err(vcmax_Tsen),sampleNo=noNAlength(vcmax_Tsen)))
write.csv(vcmax_Tsen, file='outputLeafOnto/vcmax_Tsen.csv', row.names=F)

summary(lm(vcmax_Tsen~spp*leafAge, data=jun))
anova(lm(vcmax_Tsen~spp*leafAge, data=jun))
summary(lm(vcmax_Tsen~spp*treatment, data=aug))
anova(lm(vcmax_Tsen~spp*treatment, data=aug))

vcmax_TsenShort <- as.data.frame(summarise(group_by(subset(gasex, roundNum==1), spp, leafAge, treatment), 
                                 vcmax_TsenMean=mean.na(vcmax_Tsen),vcmax_TsenSE=s.err(vcmax_Tsen),sampleNo=noNAlength(vcmax_Tsen)))
write.csv(vcmax_TsenShort, file='outputLeafOnto/vcmax_TsenShort.csv', row.names=F)

summary(lm(vcmax_Tsen~spp*leafAge, data=subset(jun, roundNum==1)))
anova(lm(vcmax_Tsen~spp*leafAge, data=subset(jun, roundNum==1)))
summary(lm(vcmax_Tsen~spp*treatment, data=subset(aug, roundNum==1)))
anova(lm(vcmax_Tsen~spp*treatment, data=subset(aug, roundNum==1)))

keep <- read.csv('dataLeafOnto/keepLocalMax.csv')
jun <- merge(jun, keep[,c('rep','keepMin','keepMax')], by='rep', all.x=T, all.y=F)
aug <- merge(aug, keep[,c('rep','keepMin','keepMax')], by='rep', all.x=T, all.y=F)
jun$time2 <- ifelse(jun$keepMax=='yes', 'midmorning', 'no')
jun$time2 <- ifelse(jun$keepMin=='yes', 'midday', jun$time2)
aug$time2 <- ifelse(aug$keepMax=='yes', 'midmorning', 'no')
aug$time2 <- ifelse(aug$keepMin=='yes', 'midday', aug$time2)

summary(lm(vcmax_Tsen~spp*leafAge*time2, data=subset(jun, time2=='midmorning' | time2=='midday')))
anova(lm(vcmax_Tsen~spp*leafAge*time2, data=subset(jun, time2=='midmorning' | time2=='midday')))
summary(lm(vcmax_Tsen~spp*treatment*time2, data=subset(aug, time2=='midmorning' | time2=='midday')))
anova(lm(vcmax_Tsen~treatment*time2*spp, data=subset(aug, time2=='midmorning' | time2=='midday')))
