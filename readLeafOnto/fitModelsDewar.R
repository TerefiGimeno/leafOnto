#function to fit the CAP version from Dewar et al. 2017
#Gamma* is 42.75, from Bernacchi et al. 2001
fitUSO <- function(df){
  myFit <- nls(gc ~ (1+Ep/sqrt(VpdL))*(Photo/(CO2S-42.75)), start=list(Ep=3), data=df)
  return(myFit)
}
fitMES <- function(df){
  myFit <- lm(gc~MESindex-1, data=df)
  return(myFit)
}
gasexL <- list()
junM <- subset(gasex, campaign=='Jun-15' & leafAge=='mature')
junY <- subset(gasex, campaign=='Jun-15' & leafAge=='young')
augMC <- subset(gasex, campaign=='Aug-15' & leafAge=='mature' & treatment=='control')
augMD <- subset(gasex, campaign=='Aug-15' & leafAge=='mature' & treatment=='drought')
#augYC <- subset(gasex, campaign=='Aug-15' & leafAge=='young' & treatment=='control')
#augYD <- subset(gasex, campaign=='Aug-15' & leafAge=='young' & treatment=='drought')
mySpp <- c(levels(gasex$spp))
for (i in 1:length(mySpp)){
  gasexL[[i]] <- subset(junM, spp==mySpp[i])
}
for (i in 1:length(mySpp)){
  gasexL[[i+length(mySpp)]] <- subset(junY, spp==mySpp[i])
}
for (i in 1:length(mySpp)){
  gasexL[[i+2*length(mySpp)]] <- subset(augMC, spp==mySpp[i])
}
for (i in 1:length(mySpp)){
  gasexL[[i+3*length(mySpp)]] <- subset(augMD, spp==mySpp[i])
}
#mySpp2 <- c('PTE','PIN','BET','EUC','MOL','ZEA')
#for (i in 1:length(mySpp2)){
# gasexL[[i+4*length(mySpp)]] <- subset(augYC, spp==mySpp2[i])
#}
#for (i in 1:length(mySpp2)){
# gasexL[[i+4*length(mySpp)+length(mySpp2)]] <- subset(augYD, spp==mySpp2[i])
#}

modelFits <- lapply(gasexL, fitUSO)
modelFitsMES <- lapply(gasexL, fitMES)
resultsFits <- data.frame(row.names = 1:length(gasexL))
for (i in 1:length(gasexL)){
  resultsFits$Campaign[i] <- gasexL[[i]][1, 'campaign']
  resultsFits$Species[i] <- gasexL[[i]][1, 'spp']
  resultsFits$leafAge[i] <- gasexL[[i]][1, 'leafAge']
  resultsFits$Treatment[i] <- gasexL[[i]][1, 'treatment']
  resultsFits$Ep[i] <- coefficients(modelFits[[i]])[1]
  resultsFits$Ep_err[i] <- summary(modelFits[[i]])$coefficients[2]
  resultsFits$t_val_Ep[i] <- summary(modelFits[[i]])$coefficients[3]
  resultsFits$P_val_Ep[i] <- summary(modelFits[[i]])$coefficients[4]
  resultsFits$upCI_Ep[i] <- confint(modelFits[[i]])[1]
  resultsFits$lwCI_Ep[i] <- confint(modelFits[[i]])[2]
  resultsFits$AIC_CAP[i] <- AIC(modelFits[[i]])
  resultsFits$RMSE_CAP[i] <- sqrt(sum((gasexL[[i]][,'gc']-predict(modelFits[[i]]))^2)
                           /noNAlength(gasexL[[i]][,'gc']))
  resultsFits$sqrtEmax[i] <- coefficients(modelFitsMES[[i]])[1]
  resultsFits$sqrtEmax_err[i] <- summary(modelFitsMES[[i]])$coefficients[2]
  resultsFits$t_val_MES[i] <- summary(modelFitsMES[[i]])$coefficients[3]
  resultsFits$P_val_MES[i] <- summary(modelFitsMES[[i]])$coefficients[4]
  resultsFits$AIC_MES[i] <- AIC(modelFitsMES[[i]])
  resultsFits$R2_MES[i] <- summary(modelFitsMES[[i]])$r.squared
  resultsFits$RMSE_MES[i] <- sqrt(sum((gasexL[[i]][which(!is.na(gasexL[[i]][,'MESindex'])),'gc']-
                                      predict(modelFitsMES[[i]]))^2)/noNAlength(gasexL[[i]][,'gc']))
}
resultsFits$Emax <- resultsFits$sqrtEmax^2
resultsFits$Emax_err <- resultsFits$sqrtEmax_err^2
resultsFits$Campaign <- ifelse(resultsFits$Campaign==1,'Aug-2015','Jun-2015')
resultsFits$leafAge <- ifelse(resultsFits$leafAge==1,'mature','young')
resultsFits$Treatment <- ifelse(resultsFits$Treatment==2,'drought','well-watered')
mySpp <- data.frame(cbind(c(levels(gasex$spp)), c(1:length(levels(gasex$spp)))))
names(mySpp) <- c('spp', 'Species')
resultsFits <- merge(resultsFits, mySpp, by='Species', all=T)
source('readLeafOnto/morpho.R')
resultsFits <- merge(resultsFits, morpho, by=c('spp','leafAge'), all=T)
#the ANOVAS and homogeneity of slopes
anova(lm(Ep~leafAge, data=subset(resultsFits, Campaign=='Jun-2015')))
summary(lm(Ep~leafAge*lmaMean, data=subset(resultsFits, Campaign=='Jun-2015')))
summary(lm(Ep~leafAge*thickMean, data=subset(resultsFits, Campaign=='Jun-2015')))
anova(lm(Emax~leafAge, data=subset(resultsFits, Campaign=='Jun-2015')))
summary(lm(Emax~leafAge*lmaMean, data=subset(resultsFits, Campaign=='Jun-2015')))
summary(lm(Emax~leafAge*thickMean, data=subset(resultsFits, Campaign=='Jun-2015')))
anova(lm(Ep~Treatment, data=subset(resultsFits, Campaign=='Aug-2015')))
anova(lm(Emax~Treatment, data=subset(resultsFits, Campaign=='Aug-2015')))

anova(lm(Ep~lmaMean, data=subset(resultsFits, Campaign=='Jun-2015')))
anova(lm(Emax~lmaMean, data=subset(resultsFits, Campaign=='Jun-2015')))
anova(lm(Ep~thickMean, data=subset(resultsFits, Campaign=='Jun-2015')))
anova(lm(Emax~thickMean, data=subset(resultsFits, Campaign=='Jun-2015')))

figureS2 <- subset(resultsFits,Campaign==
                        'Jun-2015')[,c('spp','leafAge','Ep','Ep_err',
                                       'Emax','Emax_err','lmaMean','lmaSE','thickMean','thcikSE')]
figureS2 <- reshape(figureS2, v.names=c('leafAge','Ep','Ep_err',
                                              'Emax','Emax_err','lmaMean','lmaSE','thickMean','thcikSE'), 
                       idvar='spp', timevar='leafAge', direction='wide')
write.csv(figureS2, file='outputLeafOnto/figureS2.csv', row.names=F)
