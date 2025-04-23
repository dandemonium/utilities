;+
; NAME:
;   STROM_CONV
; PURPOSE:
;    Translates stromgren color combinations from the catalog to
;    individual uvby magnitudes
; Modification 
;    2018-04-12: Jason Eastman, CfA
;                Renamed, documented, and cleaned up for distribution with EXOFASTv2
;    2023-02-22: Dan Stevens, UM-Duluth
;		  Switch from Gaia EDR3 to DR3
;-
function strom_conv,V,sigV,by,sigby,m1,sigm1,c1,sigc1,silent=silent,useticav=useticav

if n_params() lt 8 then begin
  print,'syntax: result=strom_conv(V,sigV,by,sigby,m1,sigm1,c1,sigc1)'
  retall
endif

y = V
b = V + by 
u = V + 3*by + 2*m1 + c1
v = V + 2*by + m1

sigy = sigV
sigb = sqrt(sigV^2 + sigby^2)
sigu = sqrt(sigV^2 + (3*sigby)^2 + (2*sigm1)^2 + sigc1^2)
sigv = sqrt(sigV^2 + (2*sigby)^2 + sigm1^2)

if not keyword_set(silent) then begin
  print,'uvby:',u,v,b,y
  print,'sig_uvby:',sigu,sigv,sigb,sigy
endif

return,[u,sigu,v,sigv,b,sigb,y,sigy]

end

pro mkticsed_ds, ticid, priorfile=priorfile, rvfile=rvfile, france=france, ra=ra, dec=dec

if !version.os_family eq 'Windows' then $
   message,'This program relies on queryvizier, which is not supported in Windows'

;; for use without a license
if lmgr(/vm) or lmgr(/runtime) then begin
   par = command_line_args(count=numargs)
   
   ticid = strtrim(par[0],2)
   for i=0L, numargs-1 do begin
      if strpos(par[i],'=') ne -1 then begin
         entries = strsplit(par[i],'=',/extract)
         if strupcase(entries[0]) eq 'PRIORFILE' then priorfile = strtrim(entries[1],2)
         if strupcase(entries[0]) eq 'SEDFILE' then sedfile = strtrim(entries[1],2)
         if strupcase(entries[0]) eq 'FRANCE' then france = long(entries[1])
         if strupcase(entries[0]) eq 'RA' then ra = double(entries[1])
         if strupcase(entries[0]) eq 'DEC' then dec = double(entries[1])
         if strupcase(entries[0]) eq 'RVFILE' then rvfile = strtrim(entries[1],2) 
      endif
   endfor
endif

if ~keyword_set(rvfile) then rvfile='_'+string(ticid)+'.apogee.dr17.rv'
if keyword_set(france) then cfa = 0B $
else cfa = 1B

;; match RA/Dec to TICv8 (closest)
if n_elements(ra) ne 0 and n_elements(dec) ne 0 then begin
   print, 'WARNING: querying by RA/Dec is less robust than querying by TIC ID and may lead to misidentification'   
   qtic = queryvizier('IV/38/tic',[ra,dec],2d0, /silent, /all, cfa=cfa)
   if (size(qtic))[2] ne 8 then begin
      print, 'No match to ' + string(ra,dec,format='(f0.8,",",f0.8)')
      return
   endif
   junk = min(qtic._r,ndx)
   ticid = qtic[ndx].tic
   print, 'Matching TIC ID is ' + strtrim(ticid,2)
endif

if n_elements(ticid) eq 0 then message, 'TICID is required'
if n_elements(priorfile) eq 0 then priorfile = ticid + '.priors'
if n_elements(sedfile) eq 0 then sedfile = ticid + '.sed'

dist = 120d0

;; query TICv8 for the TIC ID
if strtrim(long(ticid),2) eq ticid then begin
   qtic = queryvizier('IV/38/tic','TIC ' + strtrim(ticid,2),/allcolumns,cfa=cfa)
endif else begin
   ;; query by supplied name (less robust)
   print, 'WARNING: querying by name is less robust than querying by TIC ID and may lead to misidentification'
   qtic = queryvizier('IV/38/tic',ticid,2d0,/allcolumns,cfa=cfa)
   qgaia = queryvizier('I/345/gaia2',ticid,2d0,/allcolumns,cfa=cfa)   

   sorted = sort(qgaia.gmag)

   nstars = n_elements(qgaia)
   if strpos(ticid,'B') eq (strlen(ticid)-1) and nstars ge 2 then begin
      gaiaid = qgaia[sorted[1]].source
   endif else gaiaid = qgaia[sorted[0]].source
   match = where(qtic.gaia eq gaiaid)
   if match[0] eq -1 then message, 'no matching star found. try using the TIC ID directly'
   qtic = qtic[match]

   ticid = strtrim(qtic.tic,2)
   print, 'Matching TIC ID is ' + strtrim(ticid,2)
endelse

if (size(qtic))[2] ne 8 then begin
   print, 'No match to ' + strtrim(ticid,2)
   return
endif

match = (where(qtic.tic eq ticid))[0]
if match eq -1 then begin
   print, 'No match to ' + strtrim(ticid,2)
   return
endif

qtic = qtic[match]
star = [qtic.raj2000,qtic.dej2000]

;; prior file
openw, priorlun, priorfile, /get_lun
printf, priorlun, '#### TICv8 ####'
if finite(qtic.mass) and finite(qtic.rad) and finite(qtic.teff) then begin
   ;; require all three. starting hybrid
   ;; values is less robust than starting at solar
   printf, priorlun, qtic.mass, format='("mstar",x,f0.2)'
   printf, priorlun, qtic.rad, format='("rstar",x,f0.2)'
   printf, priorlun, qtic.teff, format='("teff",x,i5)'
endif

;; use the TICv8 prior on metallicity if available
feh = qtic._m_h_
;; this order preserves NaN, sets floor of 0.08
ufeh = 0.08d0 > qtic.e__m_h_ 
if finite(feh) and finite(ufeh) then begin
   printf, priorlun, feh, '#',ufeh, format='("feh",x,f0.5,x,f0.5)'
endif

;; extinction prior 
;; Use TICv8 Gaussian prior if availble
;; upper limit from the Schlegel dust map if not
av = qtic.e_b_v*3.1d0
uav = (0.02d0 > qtic.s_e_b_v)*3.1d0
if finite(av) and finite(uav) and keyword_set(useticav) then begin
   printf, priorlun, av, uav, format='("av",x,f0.5,x,f0.5)'
   printf, priorlun, '##############'
endif else begin
   printf, priorlun, '##############'
   junk = getavprior(ra=qtic.raj2000, dec=qtic.dej2000, line=line)
   printf, priorlun, line
endelse

;if finite(qtic.teff) then printf, priorlun, qtic.teff, format='("teff",x,i5)'

;; open the SED file for writing
;fmt = '(a13,x,f9.6,x,f0.6,x,f0.6)'
;openw, lun, sedfile, /get_lun
;printf, lun, '# bandname magnitude used_errors catalog_errors'

;; use the Gaia ID to query the Gaia catalog
gaiaid = qtic.gaia

;; DR3 (print BP/RP/G, but leave it commented out)
;qgaia3=queryvizier('I/350/gaiaedr3',star,dist/60.,/silent,cfa=cfa,/all)
qgaia3=queryvizier('I/355/gaiadr3',star,dist/60.,/silent,cfa=cfa,/all)
if (size(qgaia3))[2] eq 8 then begin
   match = (where(qgaia3.source eq gaiaid))[0]
   if match ne -1 then begin
      qgaia3 = qgaia3[match]
      
      if finite(qgaia3.plx) and finite(qgaia3.e_plx) and qgaia3.plx gt 0d0 then begin

         phot_g_mean_mag = qgaia3.gmag 
         nu_eff_used_in_astrometry = qgaia3.nueff
         pseudocolor = qgaia3.pscol
         ecl_lat = qgaia3.elat
         astrometric_params_solved = qgaia3.solved

         ;; is it within range of the Lindegren+ 2020 prescription?
         if ( (astrometric_params_solved eq 31 and nu_eff_used_in_astrometry ge 1.1d0 and nu_eff_used_in_astrometry le 1.9d0) or $
              (astrometric_params_solved eq 95 and pseudocolor ge 1.24d0 and pseudocolor le 1.72d0)) and $
            phot_g_mean_mag ge 6d0 and phot_g_mean_mag le 21d0 then begin
            zpt = get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolor, ecl_lat, astrometric_params_solved)
            printf, priorlun, "# NOTE: the Gaia DR3 parallax (" + strtrim(qgaia3.plx,2) + ") has been corrected by subtracting " + strtrim(zpt,2) + " mas according to the Lindegren+ 2020 prescription"
            printf, priorlun, qgaia3.plx-zpt, qgaia3.e_plx, format='("parallax",x,f0.5,x,f0.5)'            
         endif else begin
            printf, priorlun, "# NOTE: the Gaia DR3 parallax could not be corrected and is raw from the catalog"
            printf, priorlun, qgaia3.plx, qgaia3.e_plx, format='("#parallax",x,f0.5,x,f0.5)'
         endelse



      endif      
      if qgaia3.gmag gt -9 and finite(qgaia3.e_gmag) and (qgaia3.e_gmag lt 1d0)  then printf, lun,'#Gaia_G_DR3',qgaia3.gmag,max([0.02d,qgaia3.e_gmag]),qgaia3.e_gmag, format=fmt
      if qgaia3.bpmag gt -9 and finite(qgaia3.e_bpmag) and (qgaia3.e_bpmag lt 1d0) then printf, lun,'#Gaia_BP_DR3',qgaia3.bpmag,max([0.02d,qgaia3.e_bpmag]),qgaia3.e_bpmag, format=fmt
      if qgaia3.rpmag gt -9 and finite(qgaia3.e_rpmag) and (qgaia3.e_rpmag lt 1d0) then printf, lun,'#Gaia_RP_DR3',qgaia3.rpmag,max([0.02d,qgaia3.e_rpmag]),qgaia3.e_rpmag, format=fmt      
   endif
endif

;; use the 2MASS ID to query the 2MASS catalog
;tmassid = qtic._2mass
;q2mass=queryvizier('II/246/out',star,dist/60.,/silent,cfa=cfa)
;if (size(q2mass))[2] eq 8 then begin
;   match = (where(q2mass._2mass eq tmassid))[0]
;   if match ne -1 then begin
;      q2mass = q2mass[match]
;      if q2mass.Jmag gt -9 and (q2mass.e_Jmag lt 1d0) then printf, lun,'J2M',q2mass.Jmag,max([0.02d,q2mass.e_Jmag]),q2mass.e_Jmag, format=fmt
;      if q2mass.Hmag gt -9 and (q2mass.e_Hmag lt 1d0) then printf, lun,'H2M',q2mass.Hmag,max([0.02d,q2mass.e_Hmag]),q2mass.e_Hmag, format=fmt
;      if q2mass.Kmag gt -9 and (q2mass.e_Kmag lt 1d0) then printf, lun,'K2M',q2mass.Kmag,max([0.02d,q2mass.e_Kmag]),q2mass.e_Kmag, format=fmt
;   endif
;endif

;; use the WISE ID to query the wise catalog
;wiseid = qtic.wisea
;qwise=queryvizier('II/328/allwise',star,dist/60.,/silent,cfa=cfa)
;if (size(qwise))[2] eq 8 then begin
;   match = (where(qwise.allwise eq wiseid))[0]
;   if match ne -1 then begin
;      qwise = qwise[match]
 ;     if qwise.w1mag gt -9 and finite(qwise.e_w1mag) and (qwise.e_w1mag lt 1d0) then printf, lun,'WISE1',qwise.w1mag,max([0.03d,qwise.e_w1mag]),qwise.e_w1mag, format=fmt
 ;;     if qwise.w2mag gt -9 and finite(qwise.e_w2mag) and (qwise.e_w2mag lt 1d0) then printf, lun,'WISE2',qwise.w2mag,max([0.03d,qwise.e_w2mag]),qwise.e_w2mag, format=fmt
 ;     if qwise.w3mag gt -9 and finite(qwise.e_w3mag) and (qwise.e_w3mag lt 1d0) then printf, lun,'WISE3',qwise.w3mag,max([0.03d,qwise.e_w3mag]),qwise.e_w3mag, format=fmt
 ;     if qwise.w4mag gt -9 and finite(qwise.e_w4mag) and (qwise.e_w4mag lt 1d0) then printf, lun,'WISE4',qwise.w4mag,max([0.10d,qwise.e_w4mag]),qwise.e_w4mag, format=fmt
 ;  endif
;endif

;; use the Tycho ID to query the Stromgren catalog to get a metalicity prior
if ~finite(feh) or ~finite(ufeh) then begin
   tycid = qtic.tyc
   qpaunzen15=queryvizier('J/A+A/580/A23/catalog',star,dist/60.,/silent,/all,cfa=cfa)
   if (size(qpaunzen15))[2] eq 8 then begin
      paunzentycid = string(qpaunzen15.tyc1,qpaunzen15.tyc2,qpaunzen15.tyc3,format='(i04,"-",i05,"-",i01)')
      match = (where(paunzentycid eq tycid))[0]
      if match ne -1 then begin
         qpaunzen15 = qpaunzen15[match]
         
         ;; don't use these for the SED (?)
         if 0 then begin
            ubvymags=strom_conv(qpaunzen15.vmag,max([0.01d,qpaunzen15.e_vmag]),$
                                qpaunzen15.b_y,max([0.02d,qpaunzen15.e_b_y]),$
                                qpaunzen15.m1,max([0.02d,qpaunzen15.e_m1]),$
                                qpaunzen15.c1,max([0.02d,qpaunzen15.e_c1]),/silent)
;            printf, lun, '# Stromgren photometry, Paunzen, 2015'
;            printf, lun, '# http://adsabs.harvard.edu/abs/2015A%26A...580A..23P'
            
;            if (ubvymags(0) gt -9) and (ubvymags(1) lt 1d0) then printf, lun,'uStr',ubvymags(0),max([0.02d,ubvymags(1)]), format=fmt
;            if (ubvymags(2) gt -9) and (ubvymags(3) lt 1d0) then printf, lun,'vStr',ubvymags(2),max([0.02d,ubvymags(3)]), format=fmt
;            if (ubvymags(4) gt -9) and (ubvymags(5) lt 1d0) then printf, lun,'bStr',ubvymags(4),max([0.02d,ubvymags(5)]), format=fmt
;            if (ubvymags(6) gt -9) and (ubvymags(7) lt 1d0) then printf, lun,'yStr',ubvymags(6),max([0.02d,ubvymags(7)]), format=fmt
         endif
         
         b_y = qpaunzen15.b_y
         m1 = qpaunzen15.m1
         c1 = qpaunzen15.c1
                  
         ;; metalicity prior from Cassegrande+ 2011 (solar neighborhood)
         if b_y gt 0.23d0 and b_y lt 0.63d0 and $
            m1 gt 0.05d0 and m1 le 0.68d0 and $
            c1 gt 0.13d0 and c1 le 0.60d0 then begin
            
            ;; Cassegrande+ 2011, eq 2
            feh = 3.927d0*alog10(m1) - 14.459d0*m1^3 - 5.394d0*b_y*alog10(m1) + $
                  36.069d0*b_y*m1^3 + 3.537d0*c1*alog10(m1) - $
                  3.500d0*m1^3*c1 + 11.034d0*b_y - 22.780d0*b_y^2 + $
                  10.684d0*c1 - 6.759d0*c1^2 - 1.548d0
            ufeh = 0.10d0
            
 ;           printf, lun, '# Stromgren photometry, Paunzen, 2015'
 ;           printf, lun, '# http://adsabs.harvard.edu/abs/2015A%26A...580A..23P'
            printf, priorlun, '# Cassegrande+ 2011, eq 2'
            printf, priorlun, '#', feh, ufeh, format='("feh",x,"#",f0.5,x,f0.5)'
         endif else if b_y gt 0.43d0 and b_y lt 0.63d0 and $
            m1 gt 0.07d0 and m1 le 0.68d0 and $
            c1 gt 0.16d0 and c1 le 0.49d0 then begin
            
            ;; Cassegrande+ 2011, eq 3
            feh = -0.116d0*c1 - 1.624d0*c1^2 + 8.955d0*c1*b_y + $
                  42.008d0*b_y - 99.596d0*b_y^2 + 64.245d0*b_y^3 + $
                  8.928d0*c1*m1 + 17.275d0*m1 - 48.106d0*m1^2 + $
                  45.802d0*m1^3 - 8.467d0
            ufeh = 0.12d0
;            printf, lun, '# Stromgren photometry, Paunzen, 2015'
 ;           printf, lun, '# http://adsabs.harvard.edu/abs/2015A%26A...580A..23P' 
            printf, priorlun, '# Cassegrande+ 2011, eq 3'
            printf, priorlun,feh, ufeh, format='("feh",x,"#",f0.5,x,f0.5)'
         endif      
         
      endif
   endif
endif

print, "Querying Abdurro'uf+ (2022) for APOGEE DR17 data...'
qapo=queryvizier('III/286/allvis',star,dist/60.,/silent,cfa=cfa,/all)
if n_elements(qapo) gt 1 then begin ; multiple matches. Just take the closest.
  index = where(qapo[*]._r eq min(qapo[*]._r))
  qapo = qapo[index]
endif
if long(tag_exist(qapo,'VHelio',/quiet)) ne 0L then begin
  openw, rvlun, rvfile, /get_lun
  printf, rvlun, "#BJD	VHelio(m/s)	e_RV(m/s)"
  for i=0, n_elements(index)-1 do begin
    if qapo[i].VHelio gt -99 and not finite(qapo[i].VHelio,/nan) then begin
      ra = qapo[i].RAJ2000
      dec = qapo[i].DEJ2000
      bjd = utc2bjd(qapo[i].jd, ra, dec, earthobs='apo')
      printf, rvlun, bjd, qapo[i].VHelio*1000d0, qapo[i].e_RV*1000d0, format='(f, f, f)'
    endif
  endfor
endif else print, "No APOGEE DR17 match."

qapo2=Queryvizier('III/286/catalog',star,dist/60.,/silent,cfa=cfa,/all)
if n_elements(qapo2) gt 1 then begin ; multiple matches. Just take the closest.
  index = where(qapo2[*]._r eq min(qapo2[*]._r))
  qapo2 = qapo2[index]
endif
if long(tag_exist(qapo2,'Teff',/quiet)) ne 0L then begin
  if qapo2.Teff gt 0 and not finite(qapo2.Teff,/nan) then begin
    printf, priorlun, "### APOGEE DR17 (III/286/catalog) ###"
    printf, priorlun, 'teff', qapo2.Teff, qapo2.e_teff 
    printf, priorlun, 'feh', qapo2._FE_H_, qapo2.E__FE_H_ ; or _FE_H_SP?
  endif
endif else print, "No APOGEE DR17 match."

print, "Querying Gaia DR3 (I/357/tbosb1, I/360/binmass) RV orbital solution for initial guesses (NOT YET SUPPORTED)..."
qgaiadr3_rv1 = QueryVizier('I/357/tbosb1', star, dist/60., /silent, cfa=cfa, /all)
qgaiadr3_rv2 = Queryvizier('I/360/binmass', star, dist/60., /silent, cfa=cfa, /all)
;; we can do better than -5 < [Fe/H]_0 < +0.5
;if ~finite(feh) or ~finite(ufeh) then begin
;   printf, priorlun, '# Cassegrande+ 2011, Table 1'
;   feh = -0.06d0
;   ufeh = 0.25d0

;   printf, priorlun, '# wide Gaussian prior'
;   feh = 0d0
;;   ufeh = 1d0
;   printf, priorlun, feh, ufeh, format='("feh",x,"#",f0.5,x,f0.5)'
;endif

free_lun, priorlun
;free_lun, lun

end

