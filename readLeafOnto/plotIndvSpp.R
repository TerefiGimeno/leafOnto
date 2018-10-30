codes <- c('PTE','BET','QUE','EUC','PIN','MOL','ZEA')
myPlotCAP1 <- function(x, titulo){
  plot(subset(x, leafAge=='mature' &  treatment=='control')[,'gc']~
         subset(x, leafAge=='mature' & treatment=='control')[,'CAPindex'],col='cyan', ylim=c(0,max(x$gc)),
       xlim=c(-0.001,max(x[,'CAPindex'],na.rm=T)),ylab="CAP gc (mmol/m2/s)", xlab='', pch=19, main=titulo)
  points(subset(x, leafAge=='mature' & treatment=='drought')[,'gc']~
           subset(x, leafAge=='mature' & treatment=='drought')[,'CAPindex'],pch=19, col='red')
  points(subset(x, leafAge=='mature' & treatment=='well-watered')[,'gc']~
           subset(x, leafAge=='mature' & treatment=='well-watered')[,'CAPindex'],pch=19, col='blue')
  points(subset(x, leafAge=='young' &  treatment=='control')[,'gc']~
           subset(x, leafAge=='young' & treatment=='control')[,'CAPindex'],col='cyan',pch=17)
  points(subset(x, leafAge=='young' & treatment=='drought')[,'gc']~
           subset(x, leafAge=='young' & treatment=='drought')[,'CAPindex'],pch=17, col='red')
  points(subset(x, leafAge=='young' & treatment=='well-watered')[,'gc']~
           subset(x, leafAge=='young' & treatment=='well-watered')[,'CAPindex'],pch=17, col='blue')
}
myPlotCAP <- function(x, titulo){
  plot(subset(x, leafAge=='mature' &  treatment=='control')[,'gc']~
         subset(x, leafAge=='mature' & treatment=='control')[,'CAPindex'],col='cyan', ylim=c(0,max(x$gc)),
       xlim=c(-0.001,max(x[,'CAPindex'],na.rm=T)),ylab=" ", xlab='', pch=19, main=titulo)
  points(subset(x, leafAge=='mature' & treatment=='drought')[,'gc']~
           subset(x, leafAge=='mature' & treatment=='drought')[,'CAPindex'],pch=19, col='red')
  points(subset(x, leafAge=='mature' & treatment=='well-watered')[,'gc']~
           subset(x, leafAge=='mature' & treatment=='well-watered')[,'CAPindex'],pch=19, col='blue')
  points(subset(x, leafAge=='young' &  treatment=='control')[,'gc']~
           subset(x, leafAge=='young' & treatment=='control')[,'CAPindex'],col='cyan',pch=17)
  points(subset(x, leafAge=='young' & treatment=='drought')[,'gc']~
           subset(x, leafAge=='young' & treatment=='drought')[,'CAPindex'],pch=17, col='red')
  points(subset(x, leafAge=='young' & treatment=='well-watered')[,'gc']~
           subset(x, leafAge=='young' & treatment=='well-watered')[,'CAPindex'],pch=17, col='blue')
}
myPlotMES1 <- function(x){
  plot(subset(x, leafAge=='mature' &  treatment=='control')[,'gc']~
         subset(x, leafAge=='mature' & treatment=='control')[,'MESindex'],col='cyan', ylim=c(0,max(x$gc)),
       xlim=c(-0.001,max(x[,'MESindex'],na.rm=T)),ylab="MES gc (mmol/m2/s)", xlab='', pch=19)
  points(subset(x, leafAge=='mature' & treatment=='drought')[,'gc']~
           subset(x, leafAge=='mature' & treatment=='drought')[,'MESindex'],pch=19, col='red')
  points(subset(x, leafAge=='mature' & treatment=='well-watered')[,'gc']~
           subset(x, leafAge=='mature' & treatment=='well-watered')[,'MESindex'],pch=19, col='blue')
  points(subset(x, leafAge=='young' &  treatment=='control')[,'gc']~
           subset(x, leafAge=='young' & treatment=='control')[,'MESindex'],col='cyan',pch=17)
  points(subset(x, leafAge=='young' & treatment=='drought')[,'gc']~
           subset(x, leafAge=='young' & treatment=='drought')[,'MESindex'],pch=17, col='red')
  points(subset(x, leafAge=='young' & treatment=='well-watered')[,'gc']~
           subset(x, leafAge=='young' & treatment=='well-watered')[,'MESindex'],pch=17, col='blue')
}
myPlotMES <- function(x){
  plot(subset(x, leafAge=='mature' &  treatment=='control')[,'gc']~
         subset(x, leafAge=='mature' & treatment=='control')[,'MESindex'],col='cyan', ylim=c(0,max(x$gc)),
       xlim=c(-0.001,max(x[,'MESindex'],na.rm=T)),ylab=" ", xlab='', pch=19)
  points(subset(x, leafAge=='mature' & treatment=='drought')[,'gc']~
           subset(x, leafAge=='mature' & treatment=='drought')[,'MESindex'],pch=19, col='red')
  points(subset(x, leafAge=='mature' & treatment=='well-watered')[,'gc']~
           subset(x, leafAge=='mature' & treatment=='well-watered')[,'MESindex'],pch=19, col='blue')
  points(subset(x, leafAge=='young' &  treatment=='control')[,'gc']~
           subset(x, leafAge=='young' & treatment=='control')[,'MESindex'],col='cyan',pch=17)
  points(subset(x, leafAge=='young' & treatment=='drought')[,'gc']~
           subset(x, leafAge=='young' & treatment=='drought')[,'MESindex'],pch=17, col='red')
  points(subset(x, leafAge=='young' & treatment=='well-watered')[,'gc']~
           subset(x, leafAge=='young' & treatment=='well-watered')[,'MESindex'],pch=17, col='blue')
}

windows(12,8)
par(mfrow=c(2,7))
myPlotCAP1(subset(gasex, spp==codes[1]), codes[1])
legend('topright', pch=c(17,19), col='black', legend=c('Developing','Mature'), bty='n', cex=0.7)
for (i in 2:length(codes)){
  myPlotCAP(subset(gasex, spp==codes[i]), codes[i])
}
myPlotMES1(subset(gasex, spp==codes[1]))
legend('topleft', pch=c(15,15,15), legend=c('WW Jun','WW Aug', 'WS Aug'), col=c('blue', 'cyan', 'red'), bty='n', cex=0.7)
for (i in 2:length(codes)){
  myPlotMES(subset(gasex, spp==codes[i]))
}
