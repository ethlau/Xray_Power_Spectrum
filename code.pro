h=where(1./k_w gt t1 and 1./k_w lt t2)
;Computing the power in each channel  1,2,4
;and the cross powers with the X-ray
n = 17.0
theta = findgen(n)/(n-1.0)*360.0*!DtoR
x = 1.0*sin(theta)
y = 1.0*cos(theta)

usersym, x, y
ecf=1.47e+11
;DEVICE, RETAIN=2, DECOMPOSED=0
colours
flux1=(readfits('egs_flu_ch1.fits',hea1));[0:1131,*]
fab1=(readfits('ch1_ab_rep.fits',hea6));[0:1131,*]
t_ch1=(readfits('egs_exp_ch1.fits',hea3));[0:1131,*]
mean2=(readfits('mean_0520_flux.fits',hea3));[0:1131,*]
flux2=(readfits('egs_flu_ch2.fits',hea1));[0:1131,*]
fab2=(readfits('ch2_ab_rep.fits',hea6));[0:1131,*]
t_ch2=(readfits('egs_exp_ch2.fits',hea3));[0:1131,*]
flux4=(readfits('egs_ch4_flat_rep.fits',hea1));[0:1131,*]
fab4=(readfits('ch4_ab_rep.fits',hea6));[0:1131,*]
t_ch4=(readfits('egs_exp_ch4.fits',hea3));[0:1131,*]
img=(readfits('flux_rep.fits',h));[0:1131,*]
flux=(readfits('deltaf_0520_new.fits',hea2));[0:1131,*]
t_x=(readfits('exp_ab_0520.fits',hea4));[0:1131,*]
cmask=t_ch1/t_ch1;(readfits('mask_rep_corr.fits',hea5));[0:1131,*]
xab=(readfits('deltaf_ab_0520.fits',hea5));[0:1131,*]
print, mean(flux)
t_x=t_x
s_e=size(flux1)
nx_tot=s_e(1)
ny_tot=s_e(2)
;cmask=fltarr(nx_tot-1,ny_tot-1)
;mm=where(flux1 eq 0l , count)
;if count gt -1 then cmask(mm)=1
;if count gt -1 then cmask(mm)=0
;FIXING THE MASK AFTER PROJECTION
mm=where(cmask lt 1.0 , count)
if count gt -1 then cmask(mm)=0
cmask(*,0)=0.
cmask(*,1)=0.
writefits, 'cmask_fixex.fits', cmask, hea5
flux1=flux1*cmask
flux2=flux2*cmask
flux4=flux4*cmask
flux=flux*cmask/ecf
print, mean(flux)
mean2=mean2*cmask/ecf
xab=xab*cmask/ecf
fab1=fab1*cmask
fab2=fab2*cmask
fab4=fab4*cmask
s_e=size(flux1)
nx_tot=s_e(1)
ny_tot=s_e(2)

sec2rad=!pi/180/3600.
pixel=1.2 ; in arcsec
tpixel=pixel*sec2rad                            ;in radian
flux=flux/(pixel/3600.)^2*3282.8 ; from cts/s/pix--> cts/s/sr
xab=xab/(pixel/3600.)^2*3282.8
area=float(nx_tot)*float(ny_tot)*(pixel/3600.)^2/3282.8
;*tpixel^2
nx21=nx_tot/2+1
ny21=ny_tot/2+1

;Defining Binning in Real and Fourier space
exp=fltarr(nx_tot,ny_tot)
rx=fltarr(nx_tot,ny_tot)
n_k=1000
nx_center=nx21-1
ny_center=ny21-1
k_w=fltarr(nx_tot,ny_tot)
k_x=shift(dist(nx_tot,1),-nx21)
k_y=shift(dist(ny_tot,1),-ny21)
for i=0l,nx_tot-1 do k_w(i,*)=sqrt((k_x(i)/nx_tot)^2+(k_y(*)/ny_tot)^2)/pixel
k_0=1./(sqrt(2.)*nx_tot*pixel) ; 0.
k_f=1.01*max(k_w)
;tet_1grid=1.2+10^linspace(-1.,4,100)
tet_1grid=[1.2+10.^(.15*(findgen(21))),1200]
;tet_1grid=[8l, 1000l,1000l]
nq_1grid=n_elements(tet_1grid)
q_p_minmax=2.*!pi/tet_1grid
k_1_minmax=q_p_minmax/2./!pi
q_p_1grid=fltarr(nq_1grid-1)
for ik=0,nq_1grid-2 do q_p_1grid(ik)=0.5*(q_p_minmax(ik)+q_p_minmax(ik+1))

; Channel 1 (IR)

del_flux=flux1
hn0=where(cmask ne 0.,n1)
f_clip_1=float(n1)/nx_tot/ny_tot
del_flux(hn0)=del_flux(hn0)-mean(del_flux(hn0))
weight_1=fltarr(nx_tot,ny_tot)
weight_ab_1=weight_1
weight_1(hn0)=t_ch1(hn0)/mean(t_ch1(hn0))
weight_ab_1(hn0)=weight_1(hn0)
del_flux(hn0)=del_flux(hn0)*weight_1(hn0)
del_flux(hn0)=del_flux(hn0)-mean(del_flux(hn0))
del_flux_ab=fab1
del_flux_ab(hn0)=del_flux_ab(hn0)-mean(del_flux_ab(hn0))
del_flux_ab(hn0)=del_flux_ab(hn0)*weight_ab_1(hn0)
del_flux_ab(hn0)=del_flux_ab(hn0)-mean(del_flux_ab(hn0))
;remove axis in fourier space and compute fft
msk_fft=fltarr(nx_tot,ny_tot)
msk_fft(*,*)=1.
msk_fft(0,*) = 0
msk_fft(*,0) = 0
msk_fft = shift(msk_fft,-nx21,-ny21)
amp1=shift(fft(del_flux,-1,double=1),-nx21,-ny21)
ampab=shift(fft(del_flux_ab,-1,double=1),-nx21,-ny21)
amp1=amp1;;*msk_fft
ampab=ampab;;*msk_fft
;;;;;;;;;;;;;;;;;;;;;;;;;;
pairs1=dblarr(nq_1grid-1)
power1=dblarr(nq_1grid-1)
sig_p1=dblarr(nq_1grid-1)
powerab=dblarr(nq_1grid-1)
sig_pab=dblarr(nq_1grid-1)
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq)  )
if(n_elements(hp) gt 1) then begin
power1(iq)=mean(abs(amp1(hp))^2)*area/f_clip_1
pairs1(iq)=n_elements(hp)
powerab(iq)=mean(abs(ampab(hp))^2)*area/f_clip_1
endif
     endfor
writecol,'pairs_egs.dat',2.*!pi/q_p_1grid,pairs1

sig_p1=power1/(sqrt(0.5*pairs1))
sig_pab=powerab/(sqrt(0.5*pairs1))
;computing  P(A+B)-P(A-B) IR ch1
pcleanir1=power1-powerab
sig_plcir1=sqrt(sig_pab^2+sig_p1^2)

;now channel2
del_flux=flux2
hn0=where(cmask ne 0.,n1)
f_clip_1=float(n1)/nx_tot/ny_tot
del_flux(hn0)=del_flux(hn0)-mean(del_flux(hn0))
weight_2=fltarr(nx_tot,ny_tot)
weight_ab_2=weight_1
weight_2(hn0)=t_ch2(hn0)/mean(t_ch1(hn0))
weight_ab_2(hn0)=weight_2(hn0)
del_flux(hn0)=del_flux(hn0)*weight_2(hn0)
del_flux(hn0)=del_flux(hn0)-mean(del_flux(hn0))
del_flux_ab=fab2
del_flux_ab(hn0)=del_flux_ab(hn0)-mean(del_flux_ab(hn0))
del_flux_ab(hn0)=del_flux_ab(hn0)*weight_ab_2(hn0)
del_flux_ab(hn0)=del_flux_ab(hn0)-mean(del_flux_ab(hn0))
;remove axis in fourier space and compute fft
msk_fft=fltarr(nx_tot,ny_tot)
msk_fft(*,*)=1.
msk_fft(0,*) = 0
msk_fft(*,0) = 0
msk_fft = shift(msk_fft,-nx21,-ny21)
amp2=shift(fft(del_flux,-1,double=1),-nx21,-ny21)
ampab2=shift(fft(del_flux_ab,-1,double=1),-nx21,-ny21)
amp2=amp2;;*msk_fft
ampab2=ampab2;;*msk_fft
;;;;;;;;;;;;;;;;;;;;;;;;;;
pairs2=dblarr(nq_1grid-1)
power2=dblarr(nq_1grid-1)
sig_p2=dblarr(nq_1grid-1)
powerab2=dblarr(nq_1grid-1)
sig_pab2=dblarr(nq_1grid-1)
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq) )
if(n_elements(hp) gt 1) then begin
power2(iq)=mean(abs(amp2(hp))^2)*area/f_clip_1
pairs2(iq)=n_elements(hp)
powerab2(iq)=mean(abs(ampab2(hp))^2)*area/f_clip_1
endif
     endfor
sig_pab2=powerab2/(sqrt(0.5*pairs2))
;computing  P(A+B)-P(A-B) IR ch12
pcleanir2=power2-powerab2
sig_p2=power2/(sqrt(0.5*pairs2))

;now channel4
del_flux=flux4
hn0=where(cmask ne 0.,n1)
f_clip_1=float(n1)/nx_tot/ny_tot
del_flux(hn0)=del_flux(hn0)-mean(del_flux(hn0))
weight_4=fltarr(nx_tot,ny_tot)
weight_ab_4=weight_4
weight_4(hn0)=t_ch4(hn0)/mean(t_ch4(hn0))
weight_ab_4(hn0)=weight_4(hn0)
del_flux(hn0)=del_flux(hn0)*weight_4(hn0)
del_flux(hn0)=del_flux(hn0)-mean(del_flux(hn0))
del_flux_ab=fab4
del_flux_ab(hn0)=del_flux_ab(hn0)-mean(del_flux_ab(hn0))
del_flux_ab(hn0)=del_flux_ab(hn0)*weight_ab_4(hn0)
del_flux_ab(hn0)=del_flux_ab(hn0)-mean(del_flux_ab(hn0))
;remove axis in fourier space and compute fft
msk_fft=fltarr(nx_tot,ny_tot)
msk_fft(*,*)=1.
msk_fft(0,*) = 0
msk_fft(*,0) = 0
msk_fft = shift(msk_fft,-nx21,-ny21)
amp8=shift(fft(del_flux,-1,double=1),-nx21,-ny21)
ampab8=shift(fft(del_flux_ab,-1,double=1),-nx21,-ny21)

;;;;;;;;;;;;;;;;;;;;;;;;;;
pairs8=dblarr(nq_1grid-1)
power8=dblarr(nq_1grid-1)
sig_p8=dblarr(nq_1grid-1)
powerab8=dblarr(nq_1grid-1)
sig_pab8=dblarr(nq_1grid-1)
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq) )
if(n_elements(hp) gt 1) then begin
power8(iq)=mean(abs(amp8(hp))^2)*area/f_clip_1
pairs8(iq)=n_elements(hp)
powerab8(iq)=mean(abs(ampab8(hp))^2)*area/f_clip_1
endif
     endfor

;computing  P(A+B)-P(A-B) IR ch12
pcleanir8=power8-powerab8
sig_p8=power8/(sqrt(0.5*pairs8))
sig_pab8=powerab8/(sqrt(0.5*pairs8))

; now X-rays
dist1=flux
del_flux=flux
hn0=where(cmask ne 0.,nn)
f_clip_1=float(nn)/nx_tot/ny_tot
del_flux(hn0)=del_flux(hn0)-mean(del_flux(hn0))
weight_3=fltarr(nx_tot,ny_tot)
weight_3(hn0)=t_x(hn0)/mean(t_x(hn0))
dist1(hn0)=del_flux(hn0)
del_flux(hn0)=del_flux(hn0)*weight_3(hn0)
del_flux(hn0)=del_flux(hn0)-mean(del_flux(hn0))
del_flux_ab=xab
del_flux_ab(hn0)=del_flux_ab(hn0)-mean(del_flux_ab(hn0))
del_flux_ab(hn0)=del_flux_ab(hn0)*weight_3(hn0)
del_flux_ab(hn0)=del_flux_ab(hn0)-mean(del_flux_ab(hn0))
msk_fft=fltarr(nx_tot,ny_tot)
msk_fft(*,*)=1.
msk_fft(0,*) = 0
msk_fft(*,0) = 0
msk_fft = shift(msk_fft,-nx21,-ny21)
amp3=shift(fft(del_flux,-1,double=1),-nx21,-ny21)
ampab3=shift(fft(del_flux_ab,-1,double=1),-nx21,-ny21)
ampab3=ampab3;;*msk_fft
amp3=amp3;;*msk_fft
pairs3=dblarr(nq_1grid-1)
power3=dblarr(nq_1grid-1)
sig_p3=dblarr(nq_1grid-1)
powerab3=dblarr(nq_1grid-1)
sig_pab3=dblarr(nq_1grid-1)
for iq=0,nq_1grid-2 do begin
 hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq) )
if(n_elements(hp) gt 1) then begin
power3(iq)=mean(abs(amp3(hp))^2)*area/f_clip_1
powerab3(iq)=mean(abs(ampab3(hp))^2)*area/f_clip_1
pairs3(iq)=n_elements(hp)
endif
endfor
sig_pab3=powerab3/(sqrt(0.5*pairs3))
sig_p3=power3/(sqrt(0.5*pairs3))
;subtract ab
pclean=power3-powerab3
writefits, 'fft_x.fits',real_part(amp3)+imaginary(amp3);*area/f_clip_1

; now cross-power ch1 x ch2
pairs12=dblarr(nq_1grid-1)
power12=dblarr(nq_1grid-1)
sig_p12=dblarr(nq_1grid-1)
amp12=(real_part(amp2)*real_part(amp1)+imaginary(amp2)*imaginary(amp1))
amp12=amp12;;*msk_fft
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq) )
if(n_elements(hp) gt 1) then begin
power12(iq)=mean(amp12(hp))*area/f_clip_1
 pairs12(iq)=n_elements(hp)
endif
     endfor
sig_p12=(sqrt(power1*power2/pairs12))

; now cross-power ch1 x X-ray
pairs1x=dblarr(nq_1grid-1)
power1x=dblarr(nq_1grid-1)
sig_p1x=dblarr(nq_1grid-1)
amp1x=(real_part(amp3)*real_part(amp1)+imaginary(amp3)*imaginary(amp1) )
amp1x=amp1x;;*msk_fft
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq) )
if(n_elements(hp) gt 1) then begin
power1x(iq)=mean(amp1x(hp))*area/f_clip_1
 pairs1x(iq)=n_elements(hp)
endif
     endfor
sig_m1x=sqrt(0.5*power1*power3/(0.5*pairs1x))
print, power1x,sig_m1x

; now cross-power ch1 ab x X-ray
pairs1xa=dblarr(nq_1grid-1)
power1xa=dblarr(nq_1grid-1)
sig_p1xa=dblarr(nq_1grid-1)
amp1xa=(real_part(amp3)*real_part(ampab)+imaginary(amp3)*imaginary(ampab))
amp1xa=amp1xa;;;*msk_fft
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq)  )
if(n_elements(hp) gt 1) then begin
power1xa(iq)=mean(amp1xa(hp))*area/f_clip_1
 pairs1xa(iq)=n_elements(hp)
endif
     endfor
sig_m1xa=sqrt(0.5*powerab*power3/(0.5*pairs1xa))

; now cross-power ch2 x X-ray
pairs2x=dblarr(nq_1grid-1)
power2x=dblarr(nq_1grid-1)
sig_p2x=dblarr(nq_1grid-1)
amp2x=(real_part(amp3)*real_part(amp2)+imaginary(amp3)*imaginary(amp2))
amp2x=amp2x;;*msk_fft
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq))
if(n_elements(hp) gt 1) then begin
power2x(iq)=mean(amp2x(hp))*area/f_clip_1
 pairs2x(iq)=n_elements(hp)
endif
     endfor
sig_m2x=sqrt(0.5*power2*power3/(0.5*pairs2x))
print, power2x,sig_m2x
; now cross-power ch2 x X-ray
pairs2xu=dblarr(nq_1grid-1)
power2xu=dblarr(nq_1grid-1)
sig_m2xu=dblarr(nq_1grid-1)
amp2xu=(real_part(amp3)*real_part(amp2)+imaginary(amp3)*imaginary(amp2))
amp2xu=amp2xu
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq) )
if(n_elements(hp) gt 1) then begin
power2xu(iq)=mean(amp2xu(hp))*area/f_clip_1
 pairs2xu(iq)=n_elements(hp)
endif
     endfor
sig_m2xu=sqrt(0.5*power2*power3/(0.5*pairs2xu))
writefits, 'fft_ch2x.fits',amp2x*area/f_clip_1

; now cross-power ch2 x X-ray ab
pairs2xa=dblarr(nq_1grid-1)
power2xa=dblarr(nq_1grid-1)
sig_p2xa=dblarr(nq_1grid-1)
amp2xa=(real_part(amp3)*real_part(ampab2)+imaginary(amp3)*imaginary(ampab2))
amp2xa=amp2xa;;*msk_fft
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq)  )
if(n_elements(hp) gt 1) then begin
power2xa(iq)=mean(amp2xa(hp))*area/f_clip_1
 pairs2xa(iq)=n_elements(hp)
endif
     endfor
sig_p2xa=sqrt(0.5*powerab2*power3/(0.5*pairs2xa))
;;;;;;;;;;;;;;;;;;;;;;;;;;

; now cross-power ch4 x X-ray
pairs8x=dblarr(nq_1grid-1)
power8x=dblarr(nq_1grid-1)
sig_p8x=dblarr(nq_1grid-1)
amp8x=(real_part(amp3)*real_part(amp8)+imaginary(amp8)*imaginary(amp3) )
amp8x=amp8x;;*msk_fft
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq) )
if(n_elements(hp) gt 1) then begin
power8x(iq)=mean(amp8x(hp))*area/f_clip_1
 pairs8x(iq)=n_elements(hp)
endif
     endfor
sig_m8x=sqrt(0.5*power8*power3/(0.5*pairs8x))
print, power8x,sig_m8x

; now cross-power ch4 ab  x X-ray
pairs8xa=dblarr(nq_1grid-1)
power8xa=dblarr(nq_1grid-1)
sig_p8xa=dblarr(nq_1grid-1)
amp8xa=(real_part(amp3)*real_part(ampab8)+imaginary(amp3)*imaginary(ampab8))
amp8xa=amp8xa;;;*msk_fft
for iq=0l,nq_1grid-2 do begin
hp=where(k_w ge k_1_minmax(iq+1) and k_w lt k_1_minmax(iq)  )
if(n_elements(hp) gt 1) then begin
power8xa(iq)=mean(amp8xa(hp))*area/f_clip_1
 pairs8xa(iq)=n_elements(hp)
endif
     endfor
sig_m8xa=sqrt(0.5*powerab8*power3/(0.5*pairs8xa))
sc=(q_p_1grid/sec2rad)^2/2./!pi
writecol,'soft_ab.txt',  2*!pi/q_p_1grid, sc*power1xa, sc*sig_m1xa, sc*power2xa, sc*sig_p2xa ,  sc*power8xa, sc*sig_m8xa, fmt='(f,e,e,e,e,e,e)'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;Propagating errors PS;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sig_plcir8=sqrt(sig_pab8^2+sig_p8^2)
sig_plcir1=sqrt(sig_pab^2+sig_p1^2)
sig_plcir2=sqrt(sig_pab2^2+sig_p2^2)
sig_pcl=sqrt(sig_p3^2+sig_pab3^2)

;writecol, '/Users/cap/goods/newdataset/ch1_a.txt', q_p_1grid/sec2rad, 2*!pi/q_p_1grid, power1x, sig_m1x, pcleanir1, sig_plcir1, pclean, sig_pcl, fmt='(f,f,e,e,e,e,e,e)'

;writecol, '/Users/cap/goods/newdataset/ch2_a.txt', q_p_1grid/sec2rad, 2*!pi/q_p_1grid, power2x, sig_m2x, pcleanir2, sig_plcir2, pclean, sig_pcl, fmt='(f,f,e,e,e,e,e,e)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;Coherence;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;coherence Ch1- X-ray;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ape
!p.multi=[0]
!p.charsize=1.
!p.symsize=.1
SET_PLOT, 'PS'
DEVICE, FILE='C1X.ps', /COLOR, BITS=8,SCALE_FACTOR=2, YSIZE=13, Xsize=13,/landsc
;window, 4, xsize=500,ysize=400

cirx1=(power1x*power1x)/(pcleanir1*pclean)
a=(2*power1x/(pcleanir1*pclean)*sig_m1x)^2
b=((power1x^2/pcleanir1*alog(abs(pclean)))*sig_pcl)^2
c=((power1x^2/pclean*alog(abs(pcleanir1)))*sig_plcir1)^2
pippo=sqrt(a+b+c)
sig_c1=pippo;cirx1*sqrt(12/pairs1x)
 plot,2*!pi/q_p_1grid,cirx1,psym=sym(1),xran=[10,500],yrange=[.001,10],/xstyle,/xlog,/ylog,xtitle=TEXTOIDL('2\pi\q'),ytitle=TEXTOIDL('C_{1,X_{1}}(q)')
xyouts, 20, 0.800, 'Ch1 vs 0.5-2 keV'
errplot,2*!pi/q_p_1grid,cirx1-sig_c1,cirx1+sig_c1
DEVICE, /close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;coherence Ch2- X-ray;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=[0]

SET_PLOT, 'PS'
DEVICE, FILE='C2X.ps', /COLOR, BITS=8, XSIZE=13, ysize=13,/landscape
!p.charsize=1.
!p.symsize=.1
!p.charthick=3
!p.thick=2
;window, 3, xsize=500,ysize=400

cirx2=(power2x*power2x)/(pcleanir2*pclean)
a=(2*power2x/(pcleanir2*pclean)*sig_m2x)^2
b=((power2x^2/pcleanir2*alog(abs(pclean)))*sig_pcl)^2
c=((power2x^2/pclean*alog(abs(pcleanir2)))*sig_plcir2)^2
pippo=sqrt(a+b+c)
sig_c2=pippo
 plot,2*!pi/q_p_1grid,cirx2,psym=sym(1),xran=[10,1000],yrange=[.0001,1.],/xstyle,/xlog,/ylog,xtitle=TEXTOIDL('2\pi\q'),ytitle=TEXTOIDL('C_{2,X_{1}}(q)')
errplot,2*!pi/q_p_1grid,cirx2-sig_c2,cirx2+sig_c2
xyouts, 20, 0.800, 'Ch2 vs 0.5-2 keV'
DEVICE, /close
writecol, 'coh2_052.txt', 2*!pi/q_p_1grid,cirx2,sig_c2,fmt='(f,f,f)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;AB-PLOTS;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;coherence Ch4- X-ray;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=[0]

SET_PLOT, 'PS'
DEVICE, FILE='C4X.ps', /COLOR, BITS=8, XSIZE=13, ysize=13,/landscape
!p.charsize=1.
!p.symsize=.1
!p.charthick=3
!p.thick=2
;window, 3, xsize=500,ysize=400

cirx8=(power8x*power8x)/(pcleanir8*pclean)
a=(2*power8x/(pcleanir8*pclean)*sig_m8x)^2
b=((power8x^2/pcleanir2*alog(abs(pclean)))*sig_pcl)^2
c=((power8x^2/pclean*alog(abs(pcleanir8)))*sig_plcir8)^2
pippo=sqrt(a+b+c)
sig_c8=pippo
 plot,2*!pi/q_p_1grid,cirx8,psym=sym(1),xran=[10,500],yrange=[.001,10.],/xstyle,/xlog,/ylog,xtitle=TEXTOIDL('2\pi\q'),ytitle=TEXTOIDL('C_{2,X_{1}}(q)')
errplot,2*!pi/q_p_1grid,cirx8-sig_c8,cirx8+sig_c8
xyouts, 20, 0.800, 'Ch4 vs 0.5-2 keV'
DEVICE, /close
s2=(2*sig_p12/power12)*power12^2
writecol, 'coh4_052.txt', 2*!pi/q_p_1grid,cirx8,sig_c8,fmt='(f,f,f)'

writecol, 'irps.txt',2*!pi/q_p_1grid,(q_p_1grid/sec2rad)^2*pcleanir1/2/!pi,(q_p_1grid/sec2rad)^2*sig_plcir1/2/!pi,(q_p_1grid/sec2rad)^2*pcleanir2/2/!pi,(q_p_1grid/sec2rad)^2*sig_plcir2/2/!pi,(q_p_1grid/sec2rad)^2*power12/2/!pi,(q_p_1grid/sec2rad)^2*sig_p12/2/!pi,fmt='(f,f,f,f,f,f,f)'
SET_PLOT, 'PS'
;plot X-ray vs ch1
DEVICE, FILE='C1X_ps.ps', /COLOR, BITS=8,YSIZE=13, XSIZE=13,/landscape
!p.multi=[0];,3,1]
!p.charsize=1.
del1=(q_p_1grid/sec2rad)^2*pcleanir1/2./!pi
sig_del1=(q_p_1grid/sec2rad)^2*sig_plcir2/2./!pi
;plot X-ray PS + A+B and A-B
delc=(q_p_1grid/sec2rad)^2*pclean/2./!pi
sig_delc=(q_p_1grid/sec2rad)^2*sig_pcl/2./!pi

delAB=(q_p_1grid/sec2rad)^2*powerab8/2./!pi
;plot cross Power Ch1-Xray
del1x=(q_p_1grid/sec2rad)^2*power1x/2./!pi
sig_del1x=(q_p_1grid/sec2rad)^2*sig_m1x/2./!pi
del1xa=(q_p_1grid/sec2rad)^2*power1xa/2./!pi

plot_oo,2*!pi/q_p_1grid,del1x,psym=6,xran=[10,800],/xs,yran=[1e-12,1e-9],/ys,xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi (erg s^{-1} cm^{-2} sr^{-1} nW m^{-2} sr^{-1})'), Symsize=0.4
errplot,2*!pi/q_p_1grid,del1x-sig_del1x,del1x+sig_del1x
xyouts, 20, 5e-10, 'Ch1 vs 0.5-2 keV'
DEVICE, /close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;PLOTS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

SET_PLOT, 'PS'
;plot X-ray vs ch4
DEVICE, FILE='C4X_ps.ps', /COLOR, BITS=8, XSIZE=13, ysize=13,/landscape
!p.multi=[0];,3,1]
!p.charsize=1.
;window, 0, xsize=1000,ysize=400
del8=(q_p_1grid/sec2rad)^2*pcleanir8/2./!pi
sig_del8=(q_p_1grid/sec2rad)^2*sig_plcir8/2./!pi
;plot X-ray PS + A+B and A-B
delc=(q_p_1grid/sec2rad)^2*pclean/2./!pi
sig_delc=(q_p_1grid/sec2rad)^2*sig_pcl/2./!pi
delAB=(q_p_1grid/sec2rad)^2*powerab3/2./!pi
;plot cross Power Ch4-Xray
del8x=(q_p_1grid/sec2rad)^2*power8x/2./!pi
del8xa=(q_p_1grid/sec2rad)^2*power8xa/2./!pi
sig_del8x=(q_p_1grid/sec2rad)^2*sig_m8x/2./!pi
plot_oo,2*!pi/q_p_1grid,del8x,psym=6,xran=[10.,800],/xs,yran=[1e-12,1e-9],/ys,xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('q^2P_{4,X}(q)/2\pi  (erg s^{-1} cm^{-2} sr^{-1} nW m^{-2} sr^{-1})'), Symsize=0.4
errplot,2*!pi/q_p_1grid,del8x-sig_del8x,del8x+sig_del8x
xyouts, 20, 5e-10, 'Ch4 vs 0.5-2 keV'
DEVICE, /close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SET_PLOT, 'PS'
;plot X-ray vs ch2
DEVICE, FILE='C2X_ps.ps', /COLOR, BITS=8,YSIZE=13, XSIZE=13,/landscape
!p.multi=[0];,3,1]
!p.charsize=1.
del1=(q_p_1grid/sec2rad)^2*pcleanir2/2./!pi
sig_del1=(q_p_1grid/sec2rad)^2*sig_plcir2/2./!pi
;plot X-ray PS + A+B and A-B
delc=(q_p_1grid/sec2rad)^2*pclean/2./!pi
sig_delc=(q_p_1grid/sec2rad)^2*sig_pcl/2./!pi

delAB=(q_p_1grid/sec2rad)^2*powerab3/2./!pi
;plot cross Power Ch2-Xray
del2x=(q_p_1grid/sec2rad)^2*power2x/2./!pi
sig_del2x=(q_p_1grid/sec2rad)^2*sig_m2x/2./!pi
del2xa=(q_p_1grid/sec2rad)^2*power2xa/2./!pi
sig_del2xu=(q_p_1grid/sec2rad)^2*sig_m2xu/2./!pi

plot_oo,2*!pi/q_p_1grid,del2x,psym=6,xran=[10,800],/xs,yran=[1e-12,1e-9],/ys,xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('q^2P_{2,X}(q)/2\pi (erg s^{-1} cm^{-2} sr^{-1} nW m^{-2} sr^{-1})'), Symsize=0.4
errplot,2*!pi/q_p_1grid,del2x-sig_del2x,del2x+sig_del2x
xyouts, 20, 5e-10, 'Ch2 vs 0.5-7 keV'
DEVICE, /close
writecol,  '0520_ps.txt', 2*!pi/q_p_1grid,del1x, sig_del1x, del2x, sig_del2x,  cirx1,sig_c1,cirx2,sig_c2, del8x, sig_del8x,cirx8,sig_c8,del1xa,del2xa,del8xa,fmt='(f,e,e,e,e,f,f,f,f,f,f,f,f,f,f,f)'
writecol, '0520_ch.qdp',  2*!pi/q_p_1grid,cirx2, sig_c2
;DEVICE, RETAIN=2, DECOMPOSED=0
;window,0
!p.multi=0

readcol,'psgalaxiesc.dat', t,p, f='(f d)'
t=10^(t)
p=(t/sec2rad)^2*p/2/!pi
DEVICE, FILE='xrayf_s.ps', /COLOR, BITS=8,YSIZE=13, XSIZE=13,/landscape
hh=sig_delc
print,delc,sig_delc
;sig_delc4=0.5*sig_delc4/delc4*sqrt(delc4)
;delc4=sqrt(delc4)
plot_oo,2*!pi/q_p_1grid,delc,psym=6,xran=[10.,800.],/xs,yran=[1e-16,1e-11],/ys,xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('(q^2P_{X_{1}}(q)/2\pi)^{0.5} (erg cm^{-2} s^{-1} deg^{-1})'), Symsize=0.4
;oplot,1./t,sqrt(p), psym=0,col=2
errplot,2*!pi/q_p_1grid,delc-sig_delc,delc+sig_delc,col=1
readcol,'psagnyescr.dat', t,p2, f='(f d)'
t=10^(t)
p2=(t/sec2rad)^2*p2/2/!pi
oplot,1./t,sqrt(p2+p+(p+p2)/2)/3282.8, psym=0,col=4
ep=cirx2*delc
sep=sqrt((hh/delc)^2+(sig_c2/cirx2)^2)*ep
sep=0.5*sep/ep*sqrt(ep)/3282.8
ep=sqrt(ep)/3282.8
oplot,2*!pi/q_p_1grid,ep,psym=sym(2),col=0,symsize=0.4
errplot,2*!pi/q_p_1grid,ep-sep,ep+sep,col=0,symsize=0.4
writecol, 'xray_s.txt', 2*!pi/q_p_1grid,delc,sig_delc,t,p,fmt='(f,e,e,e,e)'
save, q_p_1grid,delc,sig_delc,filename='aegissoft.sav'
device, /close
!p.multi=0
!p.charsize=1.
!p.charthick=3
!p.thick=2
set_plot,"PS"
DEVICE, FILE='elements.ps', /COLOR, BITS=8,YSIZE=13, XSIZE=13
!p.multi=[0]
plot_oo,2*!pi/q_p_1grid,pairs3/2,/xs,/ys,psym=sym(12), yrange=[1,1000000],xrange=[2,1000],xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('N_q'),symsize=.5
oplot, 2*!pi/q_p_1grid,pairs3/2, col=0,psym=sym(5),symsize=.5
device, /close
set_plot, "PS"
DEVICE, FILE='counts.ps', /COLOR, BITS=8,YSIZE=13, XSIZE=13
ct=total(cmask)/(0.5*pairs3)*total(img)/total(cmask)
plot_oo,2*!pi/q_p_1grid,ct,/xs,/ys,psym=sym(5), yrange=[1,1000000],xrange=[2,1000],xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('X-ray counts'),symsize=0.4
ga=pairs2x/pairs2x*20.
s=sqrt(ct)
oplot, 2*!pi/q_p_1grid,ct, col=0,psym=sym(12),symsize=0.4
errplot, 2*!pi/q_p_1grid,ct+s,ct-2, col=0
oplot, [2,1000],[20,20],col=0
ct=total(img)/(0.5*pairs3)*total(pairs3)/total(pairs3)
s=sqrt(ct)
device, /close

set_plot,"PS"
DEVICE, FILE='accretionpower.ps', /COLOR, BITS=8,YSIZE=13, XSIZE=13,/landscape
;Effective power
colours
!p.multi=0
!p.charsize=1.
!p.charthick=3
!p.thick=2
ep=cirx2*(q_p_1grid/sec2rad)^2*pcleanir2/2./!pi
sep=sqrt((sig_del1/del1)^2+(sig_c2/cirx2)^2)*ep
;window,0
plot_oo,2*!pi/q_p_1grid,ep,psym=sym(1),xran=[4.,1000],/xs,yran=[1.e-5,1],/ys,xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('q^2P_2(q)/2\pi'), Symsize=2
errplot,2*!pi/q_p_1grid,ep-sep,ep+sep
oplot, 2*!pi/q_p_1grid, (q_p_1grid/sec2rad)^2*pcleanir2/2./!pi, psym=sym(3), col=0, Symsize=2
errplot,2*!pi/q_p_1grid,del1-sig_del1,del1+sig_del1,col=0
device, /close
pippo3=0l
pippo2=0l
for i=0, n_elements(q_p_1grid)-5 do begin
if(2*!pi/q_p_1grid(i) gt 20 and 2*!pi/q_p_1grid(i) le 200) then begin
pippo3=cirx2(i)/sig_c2(i)^2+pippo3
pippo2=(1/sig_c2(i)^2)+pippo2
endif
endfor
print, pippo3/pippo2,(1/pippo2)

t1=8.
t2=1000.
h=where(1./k_w gt t1 and 1./k_w lt t2)
a=mean(amp1x(h)  )   *area/f_clip_1
b= stdev(amp1x(h))*area/f_clip_1/sqrt(0.5*n_elements(h)  )

print,a/b

t1=8.
t2=1000.
h=where(1./k_w gt t1 and 1./k_w lt t2)
a=mean(amp2x(h)  )   *area/f_clip_1
b= stdev(amp2x(h))*area/f_clip_1/sqrt(0.5*n_elements(h)  )

print,a/b

t1=8.
t2=1000.
h=where(1./k_w gt t1 and 1./k_w lt t2)
a=mean(amp8x(h)  )   *area/f_clip_1
b= stdev(amp8x(h))*area/f_clip_1/sqrt(0.5*n_elements(h)  )

print,a/b,'amps'
print,  power1x[0],sig_m1x[0],power1x[0]/sig_m1x[0],'power ch1'
print, power2x[0], sig_m2x[0],power2x[0]/sig_m2x[0],'power ch2'
print, power8x[0], sig_m8x[0],power8x[0]/sig_m8x[0],'power ch4'

print,  power1xa[0],sig_m1xa[0],power1xa[0]/sig_m1xa[0],'power a-b ch1'
print, power2xa[0], sig_p2xa[0],power2xa[0]/sig_p2xa[0],'power a-b ch2'
print, power8xa[0], sig_m8xa[0],power8xa[0]/sig_m8xa[0],'power a-b ch4'
print, pclean[0],'X-ray'
print, pcleanir1[0],'IR'
print,cirx1[0],sig_c1[0],'coh1x'
print,cirx2[0],sig_c2[0],'coh2x'
print,cirx8[0],sig_c8[0],'coh4x'

;and ( 2*!pi/q_p_1grid gt 120 and 2*!pi/q_p_1grid le 400)))
;meanc2=total((cirx1/sig_c1^2)(where(2*!pi/q_p_1grid(i) gt 20 and
;2*!pi/q_p_1grid(i) le 90 and  2*!pi/q_p_1grid(i) gt 120 and
;2*!pi/q_p_1grid(i) le 400)))

set_plot,"PS"
DEVICE, /encapsulated,/in,xs=5,ys=4,/color,file='accretionpower.eps'
;Effective power
colours
!p.multi=0
!p.charsize=1.
!p.charthick=3
!p.thick=2
readcol,'tmk.dat', t,cr,xr,ir, f='(f d d d )'
restore,'hrk12_45mic_25mag.sav'
power2x=(q_p_1grid/sec2rad)^2*power2x/2./!pi
sig_m2x=(q_p_1grid/sec2rad)^2/2./!pi*sig_m2x
pclean=(q_p_1grid/sec2rad)^2*pclean/2./!pi
p2=(q_p_1grid/sec2rad)^2*pcleanir2/2./!pi
sir=(q_p_1grid/sec2rad)^2*sig_plcir2/2./!pi

h1=power2x^2/pclean
sh1=sqrt(  (2*power2x/pclean*sig_m2x)^2+(power2x^2/pclean^2*sig_pcl)^2)
sh1=0.5*h1^(-.5)*sh1
h1=sqrt(h1)
sir=0.5*p2^(-.5)*sir
p2=sqrt(p2)

plot_oo,2*!pi/q_p_1grid,h1,psym=sym(1),xran=[10.,2000],/xs,yran=[1.e-5,1],/ys,xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('(q^2P_2(q)/2\pi)^{0.5} (nW m^{-2} sr^{-1})'), Symsize=1
oplot,2*!pi/q_p_1grid,h1,psym=sym(1),col=2, Symsize=1
errplot, 2*!pi/q_p_1grid,h1+sh1,h1-sh1,col=2
oplot, 2*!pi/q_p_1grid, p2, psym=sym(5), col=4, Symsize=1
;oplot, t, ir,psym=0
;oplot, tpoq,DF2_DEF,col=0
;oplot, tpoq,DF2_HFE,col=0
;oplot, tpoq,DF2_LFE,col=0
errplot,2*!pi/q_p_1grid,p2-sir,p2+sir,col=4
device, /close

sig_pcl=(q_p_1grid/sec2rad)^2*sig_pcl/2./!pi
pcleanir2=(q_p_1grid/sec2rad)^2*pcleanir2/2./!pi
h1=power2x^2/pcleanir2
sh1=sqrt(  (2*power2x/pcleanir2*sig_m2x)^2+(power2x^2/pcleanir2^2*sig_plcir2)^2)

pclean=pclean
sig_pcl=sig_pcl
sh1=0.5*h1^(-.5)*sh1
h1=sqrt(h1)
sig_pcl=0.5*pclean^(-.5)*sig_pcl
pclean=sqrt(pclean)
readcol,'psgalaxiescu.dat', t,pp, f='(f d)'
t=10^(-t)*sec2rad
p2=sqrt((1./t)^2*pp/2/!pi)
device,/encapsulated,/in,xs=5,ys=4,/color,file='trasmx.eps',/landscape
plot_oo,2*!pi/q_p_1grid,h1,psym=sym(1),xran=[10.,2000],/xs,yran=[1e-12,1e-7],/ys,xtitle=TEXTOIDL('2\pi/q  (arcsec)'),ytitle=TEXTOIDL('(q^2P_2(q)/2\pi)^{0.5} (erg cm^{-2} s^{-1} sr^{-1})'), Symsize=1
oplot,2*!pi/q_p_1grid,h1,col=2,psym=sym(1),Symsize=1
errplot,2*!pi/q_p_1grid,h1+sh1,h1-sh1,col=2
oplot,2*!pi/q_p_1grid,pclean,psym=sym(5), col=4, Symsize=1
errplot,2*!pi/q_p_1grid,pclean+sig_pcl,pclean-sig_pcl,col=4
oplot,2*!pi*t/sec2rad,p2,psym=0,col=0, linestyle=2
;oplot, [10,2000],[1.6e-08,1.6e-08],psym=0
;oplot, [10,2000],[1.2e-08,1.2e-08],psym=0
readcol,'psagnyescr.dat', t,pp, f='(f d)'
t=10^(-t)*sec2rad
p3=sqrt((1./t)^2*pp/2/!pi)
oplot,2*!pi*t/sec2rad,p3,psym=0,col=0, linestyle=1
oplot,2*!pi*t/sec2rad,sqrt(p2^2+p3^2),psym=0,col=0
tot=sqrt(p2^2+p3^2)
model=interpol(tot,2*!pi*t/sec2rad,2*!pi/q_p_1grid)
chi=total( (h1(where(2*!pi/q_p_1grid gt 30 and   2*!pi/q_p_1grid le 1000)) -model(where(2*!pi/q_p_1grid gt 30 and   2*!pi/q_p_1grid le 1000)))/sh1(where(2*!pi/q_p_1grid gt 30 and  2*!pi/q_p_1grid le 1000 )) )^2
device,/close
writecol,'pairs.txt',2*!pi/q_p_1grid,pairs1x,fmt='(f,f)'
writecol, '/Users/cap/UDS/ele_egs.dat',2*!pi/q_p_1grid,pairs1x/2.,fmt='(f,f)'
th=2*!pi/q_p_1grid

save,filename='data_uds.xdr',th,del1x,sig_del1x,del2x,sig_del2x

end

