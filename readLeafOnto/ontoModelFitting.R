fitUSO <- function(df){
  myFit <- nls(Cond ~ g0 + 1.6*(1+g1/sqrt(VpdL))*(Photo/CO2S), start=list(g0=0, g1=3), data=df)
  return(myFit)
}
gasexL <- list()
junM <- subset(gasex, campaign=='Jun-15' & leafAge=='mature')
junY <- subset(gasex, campaign=='Jun-15' & leafAge=='young')
augMC <- subset(gasex, campaign=='Aug-15' & leafAge=='mature' & treatment=='control')
augMD <- subset(gasex, campaign=='Aug-15' & leafAge=='mature' & treatment=='drought')
augYC <- subset(gasex, campaign=='Aug-15' & leafAge=='young' & treatment=='control')
augYD <- subset(gasex, campaign=='Aug-15' & leafAge=='young' & treatment=='drought')
mySpp <- c('PTE','PIN','BET','QUE','EUC','MOL','ZEA')
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
mySpp2 <- c('PTE','PIN','BET','EUC','MOL','ZEA')
for (i in 1:length(mySpp2)){
  gasexL[[i+4*length(mySpp)]] <- subset(augYC, spp==mySpp2[i])
}
for (i in 1:length(mySpp2)){
  gasexL[[i+4*length(mySpp)+length(mySpp2)]] <- subset(augYD, spp==mySpp2[i])
}

modelFits <- lapply(gasexL, fitUSO)
resultsFits <- data.frame(row.names = 1:length(gasexL))
for (i in 1:length(gasexL)){
  resultsFits$Campaign[i] <- gasexL[[i]][1, 'campaign']
  resultsFits$Species[i] <- gasexL[[i]][1, 'spp']
  resultsFits$leafAge[i] <- gasexL[[i]][1, 'leafAge']
  resultsFits$Treatment[i] <- gasexL[[i]][1, 'treatment']
  resultsFits$g0[i] <- coefficients(modelFits[[i]])[1]
  resultsFits$g0_err[i] <- summary(modelFits[[i]])$coefficients[3]
  resultsFits$t_val_g0[i] <- summary(modelFits[[i]])$coefficients[5]
  resultsFits$P_val_g0[i] <- summary(modelFits[[i]])$coefficients[7]
  resultsFits$upCI_g0[i] <- confint(modelFits[[i]])[1]
  resultsFits$lwCI_g0[i] <- confint(modelFits[[i]])[3]
  resultsFits$g1[i] <- coefficients(modelFits[[i]])[2]
  resultsFits$g1_err[i] <- summary(modelFits[[i]])$coefficients[4]
  resultsFits$t_val_g1[i] <- summary(modelFits[[i]])$coefficients[6]
  resultsFits$P_val_g1[i] <- summary(modelFits[[i]])$coefficients[8]
  resultsFits$upCI_g1[i] <- confint(modelFits[[i]])[2]
  resultsFits$lwCI_g1[i] <- confint(modelFits[[i]])[4]
  resultsFits$AIC[i] <- AIC(modelFits[[i]])
}
resultsFits$Campaign <- ifelse(resultsFits$Campaign==1,'Aug-2015','Jun-2015')
resultsFits$leafAge <- ifelse(resultsFits$leafAge==1,'mature','young')
resultsFits$Treatment <- ifelse(resultsFits$Treatment==2,'drought','well-watered')
mySpp <- c('BET','EUC','MOL','PIN','PTE','QUE','ZEA')
for (i in 1:length(mySpp)){
  resultsFits[which(resultsFits$Species==i),'Species'] <- mySpp[i]
}
write.csv(resultsFits, file='outputLeafOnto/resultsFits.csv', row.names=FALSE)
modelFitJun <-fitUSO(jun)
jun$Residuals <- residuals(modelFitJun)
summary(aov(Residuals~spp*leafAge, data=jun))
TukeyHSD(aov(Residuals~spp*leafAge, data=jun))

modelFitAug <- fitUSO(aug)
aug$Residuals <- residuals(modelFitAug)

