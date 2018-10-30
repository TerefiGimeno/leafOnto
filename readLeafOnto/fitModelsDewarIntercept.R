#function to fit the CAP version from Dewar et al. 2017
#Gamma* is 42.75, from Bernacchi et al. 2001
fitUSO <- function(df){
  myFit <- nls(gc ~ g0+(1+xi/sqrt(Dmmol))*(Photo/(CO2S-42.75)), start=list(g0=0,xi=3), data=df)
  return(myFit)
}
fitMES <- function(df){
  myFit <- lm(gc~MESindex, data=df)
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
resultsFitsInt <- data.frame(row.names = 1:length(gasexL))
for (i in 1:length(gasexL)){
  resultsFitsInt$Campaign[i] <- gasexL[[i]][1, 'campaign']
  resultsFitsInt$Species[i] <- gasexL[[i]][1, 'spp']
  resultsFitsInt$leafAge[i] <- gasexL[[i]][1, 'leafAge']
  resultsFitsInt$Treatment[i] <- gasexL[[i]][1, 'treatment']
  resultsFitsInt$g0[i] <- coefficients(modelFits[[i]])[1]
  resultsFitsInt$g0_err[i] <- summary(modelFits[[i]])$coefficients[3]
  resultsFitsInt$t_val_g0[i] <- summary(modelFits[[i]])$coefficients[5]
  resultsFitsInt$P_val_g0[i] <- summary(modelFits[[i]])$coefficients[7]
  resultsFitsInt$xi[i] <- coefficients(modelFits[[i]])[2]
  resultsFitsInt$xi_err[i] <- summary(modelFits[[i]])$coefficients[4]
  resultsFitsInt$t_val_xi[i] <- summary(modelFits[[i]])$coefficients[6]
  resultsFitsInt$P_val_xi[i] <- summary(modelFits[[i]])$coefficients[8]
  resultsFitsInt$upCI_xi[i] <- confint(modelFits[[i]])[2]
  resultsFitsInt$lwCI_xi[i] <- confint(modelFits[[i]])[4]
  resultsFitsInt$CI_xi[i] <- resultsFitsInt[i,'lwCI_Ep']-resultsFitsInt[i,'Ep']
  resultsFitsInt$AIC_CAP[i] <- AIC(modelFits[[i]])
  resultsFitsInt$RMSE_CAP[i] <- sqrt(sum((gasexL[[i]][,'gc']-predict(modelFits[[i]]))^2)
                                  /noNAlength(gasexL[[i]][,'gc']))
  resultsFitsInt$intMES[i] <- coefficients(modelFitsMES[[i]])[1]
  resultsFitsInt$intMES_err[i] <- summary(modelFitsMES[[i]])$coefficients[3]
  resultsFitsInt$t_val_intMES[i] <- summary(modelFitsMES[[i]])$coefficients[5]
  resultsFitsInt$P_val_intMES[i] <- summary(modelFitsMES[[i]])$coefficients[7]
  resultsFitsInt$sqrtEmax[i] <- coefficients(modelFitsMES[[i]])[2]
  resultsFitsInt$sqrtEmax_err[i] <- summary(modelFitsMES[[i]])$coefficients[4]
  resultsFitsInt$t_val_MES[i] <- summary(modelFitsMES[[i]])$coefficients[6]
  resultsFitsInt$P_val_MES[i] <- summary(modelFitsMES[[i]])$coefficients[8]
  resultsFitsInt$AIC_MES[i] <- AIC(modelFitsMES[[i]])
  resultsFitsInt$R2_MES[i] <- summary(modelFitsMES[[i]])$r.squared
  resultsFitsInt$RMSE_MES[i] <- sqrt(sum((gasexL[[i]][which(!is.na(gasexL[[i]][,'MESindex'])),'gc']-
                                         predict(modelFitsMES[[i]]))^2)/noNAlength(gasexL[[i]][,'gc']))
}
resultsFitsInt$Emax <- resultsFitsInt$sqrtEmax^2
resultsFitsInt$Emax_err <- 2*resultsFitsInt$sqrtEmax*resultsFitsInt$sqrtEmax_err
resultsFitsInt$g1 <- resultsFitsInt$xi*sqrt(0.1016)/1.6
resultsFitsInt$g1_err <- resultsFitsInt$xi_err*sqrt(0.1016)/1.6
resultsFitsInt$g1_CI <- (resultsFitsInt$xi-resultsFitsInt$lwCI_xi)*sqrt(0.1016)/1.6
resultsFitsInt$campaign <- ifelse(resultsFitsInt$Campaign==1,'Aug-2015','Jun-2015')
resultsFitsInt$leafAge <- ifelse(resultsFitsInt$leafAge==1,'mature','young')
resultsFitsInt$Treatment <- ifelse(resultsFitsInt$Treatment==2,'drought','well-watered')
mySpp <- data.frame(cbind(c(levels(gasex$spp)), c(1:length(levels(gasex$spp)))))
names(mySpp) <- c('spp', 'Species')
resultsFitsInt <- merge(resultsFitsInt, mySpp, by='Species', all=T)
write.csv(resultsFitsInt, file='outputLeafOnto/resultsFitsInt.csv', row.names=F)
source('readLeafOnto/morpho.R')
resultsFitsInt <- merge(resultsFitsInt, morpho, by=c('spp','leafAge'), all=T)
source('readLeafOnto/wpAnalyses.R')
anEmax <- merge(subset(resultsFitsInt, campaign=='Aug-2015'), wpSummary2, by=c('spp','Treatment'), all=T)
figure3withInt1 <- subset(resultsFitsInt, 
                          Campaign=='Jun-2015')[,c('spp','leafAge','xi2','xi2_CI')]
figure3withInt1 <- reshape(figure3withInt1, v.names=c('leafAge','xi2','xi2_CI'),
                           idvar='spp',timevar=c('leafAge'), direction='wide')
figure3withInt2 <- subset(resultsFitsInt,
                          Campaign=='Aug-2015')[,c('spp','Treatment','xi2','xi2_CI')]
figure3withInt2 <- reshape(figure3withInt2, v.names=c('Treatment','xi2','xi2_CI'),
                           idvar='spp',timevar=c('Treatment'), direction='wide')
figure3withInt <- merge(figure3withInt1, figure3withInt2, by='spp', all=T)
write.csv(figure3withInt, file='outputLeafOnto/figure3withInt.csv', row.names=F)

#the ANOVAS and homogeneity of slopes
#analyses of homegenity of slopes
table2a <- anova(lm(gc~MESindex*leafAge*spp, data=jun))
table2b <- anova(lm(gc~MESindex*treatment*spp, data=aug))
write.csv(rbind(table2a, table2b), file='outputLeafOnto/table2.csv', row.names=T)

anova(lm(Ep~leafAge, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
summary(lm(Ep~leafAge*lmaMean, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
summary(lm(Ep~leafAge*thickMean, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
anova(lm(Emax~leafAge, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
summary(lm(Emax~leafAge*lmaMean, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
summary(lm(Emax~leafAge*thickMean, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
anova(lm(Ep~Treatment, data=subset(resultsFitsInt, Campaign=='Aug-2015')))
anova(lm(Emax~Treatment, data=subset(resultsFitsInt, Campaign=='Aug-2015')))

anova(lm(Ep~lmaMean, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
anova(lm(Emax~lmaMean, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
anova(lm(Ep~thickMean, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
anova(lm(Emax~thickMean, data=subset(resultsFitsInt, Campaign=='Jun-2015')))
figureS2int <- subset(resultsFitsInt,Campaign==
                        'Jun-2015')[,c('spp','leafAge','Ep2','Ep_err',
                                       'Emax','Emax_err','lmaMean','lmaSE','thickMean','thcikSE')]
figureS2int <- reshape(figureS2int, v.names=c('leafAge','Ep','Ep_err',
                                              'Emax','Emax_err','lmaMean','lmaSE','thickMean','thcikSE'), 
                      idvar='spp', timevar='leafAge', direction='wide')
write.csv(figureS2int, file='outputLeafOnto/figureS2int.csv', row.names=F)

#plotting predicted vs. observed for CAP and MES
for (i in 1:length(gasexL)){
  gasexL[[i]][,'gc_CAP'] <- predict(modelFits[[i]],
                                    list(Photo=gasexL[[i]][,'Photo'],
                                         Dmmol=gasexL[[i]][,'Dmmol'],CO2=gasexL[[i]][,'CO2S']))
  gasexL[[i]][,'gc_MES'] <- predict(modelFitsMES[[i]],
                                    list(MESindex=gasexL[[i]][,'MESindex']))
}
gasexPred <- as.data.frame(data.table::rbindlist(gasexL))
gasexPred$qfitCAP <- ifelse(gasexPred$campaign=='Jun-15' & gasexPred$spp=='PTE', 'no', 'yes')
gasexPred$qfitCAP <- ifelse(gasexPred$campaign=='Jun-15' & gasexPred$spp=='ZEA' & gasexPred$leafAge=='young',
                            'no', gasexPred$qfitCAP)
gasexPred$qfitCAP <- ifelse(gasexPred$campaign=='Aug-15' & gasexPred$spp=='PTE' & gasexPred$treatment=='control',
                            'no', gasexPred$qfitCAP)
gasexPred$qfitCAP <- ifelse(gasexPred$campaign=='Aug-15' & gasexPred$spp=='MOL' & gasexPred$treatment=='control',
                            'no', gasexPred$qfitCAP)
gasexPred$qfitCAP <- ifelse(gasexPred$campaign=='Aug-15' & gasexPred$spp=='ZEA' & gasexPred$treatment=='control',
                            'no', gasexPred$qfitCAP)
gasexPred$qfitMES <- ifelse(gasexPred$campaign=='Jun-15' & gasexPred$spp=='PTE', 'no', 'yes')
gasexPred$qfitMES <- ifelse(gasexPred$campaign=='Aug-15' & gasexPred$spp=='PTE' & gasexPred$treatment=='control',
                            'no', gasexPred$qfitMES)
gasexPred$qfitMES <- ifelse(gasexPred$campaign=='Aug-15' & gasexPred$spp=='MOL' & gasexPred$treatment=='control',
                            'no', gasexPred$qfitMES)
summary(lm(gc_CAP~gc*leafAge, data=subset(gasexPred, qfitCAP=='yes' & campaign=='Jun-15')))
oneCAP <- lm(gc_CAP~gc, data=subset(gasexPred, qfitCAP=='yes' & campaign=='Jun-15' & leafAge=='mature'))
twoCAP <- lm(gc_CAP~gc, data=subset(gasexPred, qfitCAP=='yes' & campaign=='Jun-15' & leafAge=='young'))
summary(lm(gc_MES~gc*leafAge, data=subset(gasexPred, qfitMES=='yes' & campaign=='Jun-15')))
oneMES <- lm(gc_MES~gc, data=subset(gasexPred, qfitMES=='yes' & campaign=='Jun-15'))
df <- data.frame(row.names=1:length(oneCAP$residuals))
df$res <- abs(oneCAP$residuals)
df$model <- 'CAP'
k <- data.frame(row.names=1:length(twoCAP$residuals))
k$res <- abs(twoCAP$residuals)
k$model <- 'CAP'
df <- rbind(df, k)
k <- data.frame(row.names=1:length(oneMES$residuals))
k$res <- abs(oneMES$residuals)
k$model <- 'MES'
df <- rbind(df, k)
summary(lm(res~model, data=df))

summary(lm(gc_CAP~gc*treatment, data=subset(gasexPred, qfitCAP=='yes' & campaign=='Aug-15')))
summary(lm(gc_MES~gc*treatment, data=subset(gasexPred, qfitMES=='yes' & campaign=='Aug-15')))
threeCAP <- lm(gc_CAP~gc, data=subset(gasexPred, qfitCAP=='yes' & campaign=='Aug-15'))
threeMES <- lm(gc_MES~gc, data=subset(gasexPred, qfitMES=='yes' & campaign=='Aug-15'))
df <- data.frame(row.names=1:length(threeCAP$residuals))
df$res <- abs(threeCAP$residuals)
df$model <- 'CAP'
k <- data.frame(row.names=1:length(threeMES$residuals))
k$res <- abs(threeMES$residuals)
k$model <- 'MES'
df <- rbind(df, k)
summary(lm(res~model, data=df))


write.csv(gasexPred, file='outputLeafOnto/predVSobs.csv', row.names = F)
