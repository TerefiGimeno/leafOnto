source('readLeafOnto/fitModelsDewar.R')
names(anEmax)[2] <-'treatment'
anEmax$treatment <- ifelse(anEmax$treatment=='well-watered', 'control', anEmax$treatment)
anEmax <- merge(anEmax, gasexShortSumm[,c('spp','treatment','vcmax_TsenMean','vcmax_TsenSE')],
                by=c('spp','treatment'), all.x=T, all.y=F)
anEmax$xiSqrtVcmax0 <- anEmax$xi*sqrt(anEmax$vcmax_TsenMean)
anEmax$xiSqrtVcmax0SE <- anEmax$xi_err*sqrt(anEmax$vcmax_TsenMean)+
  0.5*anEmax$vcmax_TsenSE*anEmax$xi/sqrt(anEmax$vcmax_TsenMean)
write.csv(anEmax[,c('spp','treatment','xiSqrtVcmax0','xiSqrtVcmax0SE','xi','xi_err','wpMean','wpSE','Emax','Emax_err')],
          file='outputLeafOnto/newFigure4First.csv', row.names=F)

#plot for Belinda
anEmax$invSqrtVcmax <- sqrt(1/anEmax$vcmax_TsenMean)
anEmax$invSqrtVcmaxErr <- sqrt(1/anEmax$vcmax_TsenSE)
windows(8,8)
par(mar=c(5,6,1,1))
anEmax1 <- subset(anEmax, treatment=='control')
anEmax2 <- subset(anEmax, treatment=='drought')
Hmisc::errbar(anEmax1$invSqrtVcmax, anEmax1$xi, anEmax1$xi+anEmax1$xi_err, anEmax1$xi-anEmax1$xi_err,
              pch=19, col=as.factor(anEmax1$spp), ylab=expression(xi~(mmol^00.5~mol^-0.5)),
              xlab=expression(sqrt(1/italic(V)[cmax0])), cex=1.2, xlim=c(0.13, 0.265), ylim=c(0,11),
              cex.lab=1.3)
Hmisc::errbar(anEmax2$invSqrtVcmax, anEmax2$xi, anEmax2$xi+anEmax2$xi_err, anEmax2$xi-anEmax2$xi_err,
              pch=17, col=as.factor(anEmax1$spp), cex=1.2, add=T)
legend('bottomleft', legend=c(levels(anEmax1$spp)), pch=15, bty='n', col=1:6)
legend('topleft', legend=c('control','drought'), pch=c(19, 17), bty='n', col=1)