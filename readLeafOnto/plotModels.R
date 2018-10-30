windows(12,8)
par(mfrow=c(2,7))
mySpp <- c(levels(gasex$spp))
k <- subset(jun, leafAge=='mature' & spp==mySpp[1])
plot(k$gc~k$CAPindex, ylab='gc (mmol/m2/s)', main=mySpp[1], xlab='CAP',
     ylim=c(0,max(gasex$gc)),xlim=c(-0.01, max(gasex$CAPindex)),
     pch=19, col='blue')
p <- subset(jun, leafAge=='young' & spp==mySpp[1])
points(p$gc~p$CAPindex, pch=19, col='cyan')
r <- subset(aug, treatment=='control' & spp==mySpp[1])
points(r$gc~r$CAPindex, pch=17, col='blue')
s <- subset(aug, treatment=='drought' & spp==mySpp[1])
points(s$gc~s$CAPindex, pch=17, col='red')
for (i in 2:length(mySpp)){
  k <- subset(jun, leafAge=='mature' & spp==mySpp[i])
  plot(k$gc~k$CAPindex, ylab='', main=mySpp[i], xlab='CAP', ylim=c(0,max(gasex$gc)),
       xlim=c(-0.01, max(gasex$CAPindex)),
       pch=19, col='blue')
  p <- subset(jun, leafAge=='young' & spp==mySpp[i])
  points(p$gc~p$CAPindex, pch=19, col='cyan')
  r <- subset(aug, treatment=='control' & spp==mySpp[i])
  points(r$gc~r$CAPindex, pch=17, col='blue')
  s <- subset(aug, treatment=='drought' & spp==mySpp[i])
  points(s$gc~s$CAPindex, pch=17, col='red')
}
k <- subset(jun, leafAge=='mature' & spp==mySpp[1])
plot(k$gc~k$MESindex, ylab='gc (mmol/m2/s)', xlab='MES',
     ylim=c(0,max(gasex$gc)),xlim=c(-0.01, max(gasex$MESindex, na.rm=T)),
     pch=19, col='blue')
p <- subset(jun, leafAge=='young' & spp==mySpp[1])
points(p$gc~p$MESindex, pch=19, col='cyan')
r <- subset(aug, treatment=='control' & spp==mySpp[1])
points(r$gc~r$MESindex, pch=17, col='blue')
s <- subset(aug, treatment=='drought' & spp==mySpp[1])
points(s$gc~s$MESindex, pch=17, col='red')
legend('topleft', pch=c(19,19,17,17), col=c('blue','cyan','blue','red'),
       legend=c('WW Jun M','WW Jun D', 'WW Aug', 'WS Aug'), bty='n')
for (i in 2:length(mySpp)){
  k <- subset(jun, leafAge=='mature' & spp==mySpp[i])
  plot(k$gc~k$MESindex, ylab='', xlab='MES', ylim=c(0,max(gasex$gc)),
       xlim=c(-0.01, max(gasex$MESindex, na.rm=T)),
       pch=19, col='blue')
  p <- subset(jun, leafAge=='young' & spp==mySpp[i])
  points(p$gc~p$MESindex, pch=19, col='cyan')
  r <- subset(aug, treatment=='control' & spp==mySpp[i])
  points(r$gc~r$MESindex, pch=17, col='blue')
  s <- subset(aug, treatment=='drought' & spp==mySpp[i])
  points(s$gc~s$MESindex, pch=17, col='red')
}

