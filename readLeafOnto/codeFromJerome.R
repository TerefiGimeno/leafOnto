PRO main

@ogee_general_header.inc
s_ksi = '!9x!6'

;;... soil retention curve
b = 2.42                       ;--- unitless (Ogée & Brunet 2002, after Myrold et al. 1981)
BD = 0.25                      ;--- g/cm3 (bulk density, 4:2:1 of bark:peat:sand)
PsiSat = -35.3*9800e-6/BD^(-b) ;--- MPa (Ogée & Brunet 2002, after Myrold et al. 1981)
;;... soil/rhizosphere conductivity
KsrSat = 1e4                  ;--- mol/m2/s/MPa (Dewar et al. 2018)
;;... plant hydraulic parameters
PsiC = -2.                    ;--- MPa (Dewar et al. 2018)
Krl = [2,7,50]*1e-3           ;--- mol/m2/s/MPa (Dewar et al. 2018)
;;... plant photosynthetic parameters
Km = 710.e-6                  ;--- mol/mol (see SI)
GStar = 40e-6                 ;--- mol/mol (see SI)
Vcmax = 30.e-6                ;--- mol/m2/s (see SI)

;;... predawn range
nPsi = 100
PsiPD = -1.5*(findgen(nPsi)+0.5)/nPsi
;;...
Ksr = KsrSat*(PsiSat/PsiPD)^(2.+3./b)

;;... ksi and Emax
nSens = n_elements(Krl)
ksi  = fltarr(nPsi, nSens)
Emax = fltarr(nPsi, nSens)
for iSens=0,nSens-1 do begin
Ksl = 1./(1./Ksr + 1./Krl[iSens])
ksi[*,iSens] = sqrt(Ksl*abs(PsiC)/1.6*(Km+GStar)/Vcmax)
Emax[*,iSens] = Ksl*(PsiPD - PsiC)
endfor


;;--- open output graphic file
file = './Dewar_plot.ps'
mc_opendev,dev='psc',file=file, landscape=1

;;--- plot
!p.multi=[0,2,1]
plot,[0,0],$
  xrange=[-1.5,-0.1],xstyle=1,xtitle=s_psi+stdfont+'!dpd!n (MPa)',xminor=1,$
  yrange=[-0.4,2.2],ystyle=1,ytitle=s_ksi+stdfont+' (kPa!u0.5!n)',yminor=1,$
  /nodata,font=0,charsize=1
for iSens=0,nSens-1 do oplot,PsiPD,ksi[*,iSens], thick=5,linestyle=iSens
mc_llegend,$
  x=!x.crange[0]+(!x.crange[1]-!x.crange[0])*0.12,$
  y=!y.crange[0]+(!y.crange[1]-!y.crange[0])*0.90,$
  larr=[indgen(nSens)],$
  ltarr=[REPLICATE(5,nSens)],$
  tarr=stdfont+'K!drl!n = '+autostring(fix(Krl*1e3))+' mmol m!u-2!n s!u-1!n MPa!u-1!n',$
  row=nSens,col=1,$
  llength=3,$
  font=!p.font,/data
;;...
plot,[0,0],$
  xrange=[-1.5,-0.1],xstyle=1,xtitle=s_psi+stdfont+'!dpd!n (MPa)',xminor=1,$
  yrange=[-0.1,13],ystyle=1,ytitle=stdfont+'E!dmax!n (mmol m!u-2!n s!u-1!n)',yminor=1,$
  /nodata,font=0,charsize=1
for iSens=0,nSens-1 do oplot,PsiPD,Emax[*,iSens]*1e3, thick=5,linestyle=iSens

mc_closedev

;;--- converts to PDF and deletes PS file (for Mac only)
IF (os NE 'Windows') THEN BEGIN
cmd = 'pstopdf '+file+' -o '+(STR_SEP(file,'.ps'))[0]+'.pdf'
spawn,cmd
FILE_DELETE, file
ENDIF

stop 
END
