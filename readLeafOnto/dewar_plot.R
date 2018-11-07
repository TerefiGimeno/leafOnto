#... soil retention curve
b = 2.42                       #--- unitless (Ogée & Brunet 2002, after Myrold et al. 1981)
BD = 0.25                      #--- g/cm3 (bulk density, 4:2:1 of bark:peat:sand)
#BD = 0.43                     #bulk density in g/cm3 measured in our pots in July 2015
PsiSat = -35.3*9800e-6/BD^(-b) #--- MPa (Ogée & Brunet 2002, after Myrold et al. 1981)
#... soil/rhizosphere conductivity
KsrSat = 1e4                  #--- mol/m2/s/MPa (Dewar et al. 2018)
#... plant hydraulic parameters
PsiC = -2.                    #--- MPa (Dewar et al. 2018)
Krl = c(2,7,12,30,50)*1e-3          #--- mol/m2/s/MPa (Dewar et al. 2018)
#... plant photosynthetic parameters
Km = 710.e-6                  #--- mol/mol (see SI)
GStar = 40e-6                 #--- mol/mol (see SI)
Vcmax = 1.e-6                #--- just for conversion mol to umol

#... predawn range
nPsi = 150
PsiPD = -2*(seq(nPsi)+0.5)/nPsi
#...
Ksr = KsrSat*(PsiSat/PsiPD)^(2.+3./b) #eq 5c in Dewar 2018

#... ksi and Emax
nSens =length(Krl)
ksi  = matrix(data=0,nPsi, nSens)
Emax = matrix(data=0,nPsi, nSens)
for (iSens in 1:nSens) {
  Ksl = 1./(1./Ksr + 1./Krl[iSens])
  ksi[,iSens] = sqrt(Ksl*abs(PsiC)/1.6*(Km+GStar)/Vcmax * 101.3)
  Emax[,iSens] = Ksl*(PsiPD - PsiC) * 1e3
}

ksiExport <- data.frame(row.names=1:nPsi)
ksiExport$nPsi <- PsiPD
ksiExport$predksi2 <- ksi[,1]
ksiExport$predksi7 <- ksi[,2]
ksiExport$predksi12 <- ksi[,3]
ksiExport$predksi30 <- ksi[,4]
ksiExport$predksi50 <- ksi[,5]
ksiExport$predEmax2 <- Emax[,1]
ksiExport$predEmax7 <- Emax[,2]
ksiExport$predEmax12 <- Emax[,3]
ksiExport$predEmax30 <- Emax[,4]
ksiExport$predEmax50 <- Emax[,5]

write.csv(ksiExport, file='outputLeafOnto/newFigure4second.csv', row.names=F)

  
#--- open output graphic file
pdf("./Dewar_plot.pdf",paper='A4r')
color = c("Cyan","Red","Green")
ltype = c(1,3,2)
layout(matrix(seq(2),1,2,byrow=F))
op <- par(matrix(seq(2),1,2,byrow=F),
          mar = c(4,5,2,0) + 0.1) # (bottom,left,top,right)
plot(c(0,0),
       xlab=expression(paste(psi[pd]," (MPa)")),xlim=c(-1.5,-0.1),
       ylab=expression(paste(xi," (",kPa^0.5,")")),ylim=c(-0.4,2.2))
for (iSens in 1:nSens) {
    lines(PsiPD,ksi[,iSens],lty=ltype[iSens],lwd=3)
}
#...
plot(c(0,0),
       xlab=expression(paste(psi[pd]," (MPa)")),xlim=c(-1.5,-0.1),
       ylab=expression(paste(E[max]," (mmol ",m^-2," ",s^-1,")")),ylim=c(-0.1,13))
for (iSens in 1:nSens) {
    lines(PsiPD,Emax[,iSens]*1e3,lty=ltype[iSens],lwd=3)
}
dev.off()
par(op)

