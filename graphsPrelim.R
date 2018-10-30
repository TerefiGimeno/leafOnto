windows(12,9)
par(mfrow=c(1,3))
plot(jun$Cond~jun$VpdL, pch=19, col=jun$spp, ylab=expression(italic(g)[s]~(mmol~H[2]~O~m^-2)),
     xlab=expression(italic(D)~(kPa), cex.lab=1.5))
myPal <- c('black','red','green','blue','cyan','yellow','magenta')
palette(mypal)
legend('topright', pch=19, levels(jun$spp), col=myPal)

windows(12,9)
par(mfrow=c(3,7), mar=c(4,5,1,1))
palette(c('blue','cyan'))
k <- subset(gasex, treatment=='well-watered' & spp=='PTE')
plot(k$Cond~k$USOindex, col=k$leafAge, ylab=expression(italic(g)[s]~(mmol~H[2]~O~m^-2)),xlab='', cex.lab=1.5,
     pch=19, ylim=c(0,0.4), xlim=c(0,0.12))
legend('topright', pch=19, levels(k$leafAge), col=c('blue','cyan'))

 