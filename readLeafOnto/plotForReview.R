b = 2.42                       #--- unitless (Ogée & Brunet 2002, after Myrold et al. 1981)
BD = 0.25                      #--- g/cm3 (bulk density, 4:2:1 of bark:peat:sand)
#BD = 0.43                     #bulk density in g/cm3 measured in our pots in July 2015
PsiSat = -35.3*9800e-6/BD^(-b) #--- MPa (Ogée & Brunet 2002, after Myrold et al. 1981)
#... soil/rhizosphere conductivity
KsrSat = 1e4                  #--- mol/m2/s/MPa (Dewar et al. 2018)
nPsi = 150
PsiPD = -2*(seq(nPsi)+0.5)/nPsi
#...
Ksr = KsrSat*(PsiSat/PsiPD)^(2.+3./b) #eq 5c in Dewar 2018

myKsl <- c(0.1,1,2,3)*0.01
nSens =length(myKsl)
myKrlM = matrix(data=0,nPsi, nSens)
for (iSens in 1:nSens) {
  myKrlM[,iSens] = (1./(1./myKsl[iSens] - 1./Ksr))*1000
}