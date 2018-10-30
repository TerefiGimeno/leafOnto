#... soil retention curve
b = 2.42                       #--- unitless (Ogée & Brunet 2002, after Myrold et al. 1981)
BD = 0.25                      #--- g/cm3 (bulk density, 4:2:1 of bark:peat:sand)
PsiSat = -35.3*9800e-6/BD^(-b) #--- MPa (Ogée & Brunet 2002, after Myrold et al. 1981)
#... soil/rhizosphere conductivity
KsrSat = 1e4                  #--- mol/m2/s/MPa (Dewar et al. 2018)
#... plant hydraulic parameters
PsiC = -2.                    #--- MPa (Dewar et al. 2018)
Krl = c(2,7,50)*1e-3          #--- mol/m2/s/MPa (Dewar et al. 2018)
#... plant photosynthetic parameters
Km = 710.e-6                  #--- mol/mol (see SI)
GStar = 40e-6                 #--- mol/mol (see SI)
Vcmax = 30.e-6                #--- mol/m2/s (see SI)

#... predawn range
nPsi = 100
PsiPD = -1.5*(seq(nPsi)+0.5)/nPsi
#...
Ksr = KsrSat*(PsiSat/PsiPD)^(2.+3./b)

#... ksi and Emax
nSens =length(Krl)
ksi  = matrix(data=0,nPsi, nSens)
Emax = matrix(data=0,nPsi, nSens)
for (iSens in 1:nSens) {
    Ksl = 1./(1./Ksr + 1./Krl[iSens])
    ksi[,iSens] = sqrt(Ksl*abs(PsiC)/1.6*(Km+GStar)/Vcmax)
    Emax[,iSens] = Ksl*(PsiPD - PsiC)
}
  
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

