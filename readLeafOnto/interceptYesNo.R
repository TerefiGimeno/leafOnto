k <- subset(aug, spp=='EUC')
kww <- subset(k, treatment=='control')
kws <- subset(k, treatment=='drought')
modelMES1euc <- lm(gc~MESindex*treatment, data=k)
modelMES2euc <- lm(gc~treatment*MESindex-1, data=k)
modelCAP1euc <- lm(gc~treatment*CAPindex, data=k)
modelCAP2euc <- lm(gc~treatment*CAPindex-1, data=k)

q <- subset(aug, spp=='QUE')
qww <- subset(q, treatment=='control')
qws <- subset(q, treatment=='drought')
modelMES1que <- lm(gc~MESindex*treatment, data=q)
modelMES2que <- lm(gc~treatment*MESindex-1, data=q)
modelCAP1que <- lm(gc~treatment*CAPindex, data=q)
modelCAP2que <- lm(gc~treatment*CAPindex-1, data=q)

m <- subset(jun, spp=='MOL')
mm <- subset(m, leafAge=='mature')
my <- subset(m, leafAge=='young')
modelMES1mol <- lm(gc~MESindex*leafAge, data=m)
modelMES2mol <- lm(gc~leafAge*MESindex-1, data=m)
modelCAP1mol <- lm(gc~leafAge*CAPindex, data=m)
modelCAP2mol <- lm(gc~leafAge*CAPindex-1, data=m)

windows(12,8)
par(mfrow=c(2,3))
plot(kww$gc~kww$MESindex, pch=17, col='blue', ylab='gc (mmol/m2/s)',
     xlab='MES', ylim=c(0, max(k$gc)), xlim=c(0, max(k$MESindex, na.rm=T)))
legend('topleft', pch=c(17,17), col=c('blue', 'red'), legend=c('WW','WS'), bty='n')
points(kws$gc~kws$MESindex, pch=17, col='red')
abline(lm(k$gc~k$MESindex))
abline(lm(k$gc~k$MESindex-1), lty=5)

plot(qww$gc~qww$MESindex, pch=17, col='blue', ylab=' ',
     xlab='MES', ylim=c(0, max(q$gc)), xlim=c(0, max(q$MESindex, na.rm=T)))
legend('topleft', pch=c(17,17), col=c('blue', 'red'), legend=c('WW','WS'), bty='n')
points(qws$gc~qws$MESindex, pch=17, col='red')
abline(lm(qww$gc~qww$MESindex), col='blue')
abline(lm(qws$gc~qws$MESindex), col='red')
abline(lm(qww$gc~qww$MESindex-1), col='blue', lty=5)
abline(lm(qws$gc~qws$MESindex-1), col='red', lty=5)

plot(mm$gc~mm$MESindex, pch=19, col='blue', ylab='gc (mmol/m2/s)',
     xlab='MES', ylim=c(0, max(m$gc)), xlim=c(0, max(m$MESindex, na.rm=T)))
legend('topleft', pch=c(19,19), col=c('blue', 'cyan'), legend=c('mature','young'), bty='n')
points(my$gc~my$MESindex, pch=19, col='cyan')
abline(lm(m$gc~m$MESindex))
abline(lm(m$gc~m$MESindex-1), lty=5)

plot(kww$gc~kww$CAPindex, pch=17, col='blue', ylab='gc (mmol/m2/s)',
     xlab='CAP', ylim=c(0, max(k$gc)), xlim=c(0, max(k$CAPindex, na.rm=T)))
points(kws$gc~kws$CAPindex, pch=17, col='red')
abline(lm(k$gc~k$CAPindex))
abline(lm(k$gc~k$CAPindex-1), lty=5)
legend('topleft', legend='EUC Aug (g0_CAP>0)', bty='n')

plot(qww$gc~qww$CAPindex, pch=17, col='blue', ylab='gc (mmol/m2/s)',
     xlab='CAP', ylim=c(0, max(q$gc)), xlim=c(0, max(q$CAPindex, na.rm=T)))
points(qws$gc~qws$CAPindex, pch=17, col='red')
abline(lm(qww$gc~qww$CAPindex), col='blue')
abline(lm(qws$gc~qws$CAPindex), col='red')
abline(lm(qww$gc~qww$CAPindex-1), col='blue', lty=5)
abline(lm(qws$gc~qws$CAPindex-1), col='red', lty=5)
legend('topleft', legend='QUE Aug (g0_MES>0 in WW)', bty='n')

plot(mm$gc~mm$CAPindex, pch=19, col='blue', ylab='gc (mmol/m2/s)',
     xlab='CAP', ylim=c(0, max(m$gc)), xlim=c(0, max(m$CAPindex, na.rm=T)))
points(my$gc~my$CAPindex, pch=19, col='cyan')
abline(lm(m$gc~m$CAPindex))
abline(lm(m$gc~m$CAPindex-1), lty=5)
legend('topleft', legend='MOL Jun (g0_MES>0 in mat)', bty='n')
