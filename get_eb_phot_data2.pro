pro	get_eb_phot_data2, star, sedfile=sedfile, gaiaspfile=gaiaspfile, rvfile=rvfile, france=france, ra=ra,$
      dec=dec, dist=dist, gaiaonly=gaiaonly
;;;;;;;;; README ;;;;;;;;;;;;
; Queries these catalogs    ;
; for photometry:           ;
; GALEX (with instr. offset ;
;       and aperture corr.) ;
; Tycho-2 (w/ Mann & von    ;
;          Braun filter def'n;
; Gaia DR3 (spectrophot)    ;
;2MASS (but not 6x)	     ;
; WISE			     ;
;			    ;
; Also query for APOGEE RVs ;
; Last updated: 2025-04-23  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_params() lt 1 then begin
  print,'Syntax: get_eb_phot_data,star,sedfile=sedfile,gaiaspfile=gaiaspfile, priorfile=priorfile,france=france, ra=ra, dec=dec, dist=3.0'
  print,'"star" is any VizieR-resolvable catalog name/ID and "dist" is in arcseconds (default: 3 arcsec).'
  print,'NOTE: If setting RA and Dec (in decimal format), then "star" will be ignored.'
  retall
endif

if not(keyword_set(dist)) then dist = 3.0 ; arcsec
if n_elements(star) eq 0 then message, 'star name is required'
;if n_elements(priorfile) eq 0 then priorfile = star + '.priors'
if n_elements(sedfile) eq 0 then sedfile = star + '.sed'
if n_elements(gaiaspfile) eq 0 then gaiaspfile = star + '.gaia.sed'
if keyword_set(france) then cfa = 0B else cfa = 1B

if keyword_set(sedfile) then openw, sedlun, sedfile, /get_lun
; query TIC v8.2 by RA and Dec (if specified) or star name)
if n_elements(ra) ne 0 and n_elements(dec) ne 0 then begin
   print, 'WARNING: querying by RA/Dec is less robust than querying by catalog name (TIC ID) and may lead to misidentification'   
;   qtic = exofast_queryvizier('IV/38/tic',[ra,dec],2d0, /silent, /all, cfa=cfa)
   qtic = exofast_queryvizier('IV/39/tic82',[ra,dec],2d0, /silent, /all, cfa=cfa)
   if (size(qtic))[2] ne 8 then begin
      print, 'No match to ' + string(ra,dec,format='(f0.8,",",f0.8)')
      return
   endif
   junk = min(qtic._r,ndx)
   ticid = qtic[ndx].tic
   star = ticid
   print, 'Matching TIC ID is ' + strtrim(ticid,2)
endif else begin
   print, 'Querying TIC v8.2 (Paegert et al. 2021) for TIC ID and TESS Tmag...'
   qtic=QueryVizier('IV/39/tic82',star,dist/60.,/silent,/all,cfa=cfa)
   if (size(qtic))[2] eq 8 then begin
      match = min(qtic._r,ndx)
      ticid = qtic[ndx].tic
      if match ne -1 then begin
         qtic = qtic[match]
         if qtic.tmag gt -9 and finite(qtic.tmag) and (qtic.e_tmag lt 1d0) then begin
           printf, sedlun, '# TIC v8.2 Paegert et al. 2021; IV/39/tic82)' 
           printf, sedlun,'# TESS_TESS.Red',qtic.tmag,qtic.e_tmag,qtic.e_tmag
         endif
      endif
   endif
endelse
if ~keyword_set(rvfile) then rvfile='_'+string(ticid)+'.apogee.dr17.rv'
;; Gaia spectrophotometry
if keyword_set(gaiaspfile) then begin
   qgaiasp=Exofast_Queryvizier('I/355/xpsample',star,dist/60.,/silent,cfa=cfa,/all)
   if (size(qgaiasp))[2] eq 8 then begin
         match = min(qgaiasp._r,ndx)
      if match[0] ne -1 then begin
         source = qgaiasp[ndx].source
         match = where(qgaiasp.source eq source)
        ; if qgaiasp[match].flux[0] gt 0d0 and not finite(qgaiasp[match].e_flux[0],/nan) then begin
            ;; Gaia lambda in nm, Gaia flux in W/m^2/Hz
            ;; SED plot in microns, erg/s/cm^2
            ;; interpolate onto atmosphere wavelength scale now?
         flux = qgaiasp[match].flux*1d3*qgaiasp[match].lambda ;; W/m^2/Hz -> erg/s/cm^2
         fluxerr = qgaiasp[match].e_flux*1d3*qgaiasp[match].lambda ;; W/m^2/Hz -> erg/s/cm^2
         lambda = qgaiasp[match].lambda/1d3 ;; nm -> um
         exofast_forprint, lambda, flux, fluxerr, textout=gaiaspfile,/nocomment,format='(f5.3,x,e12.6,x,e12.6)'
;         endif else print, "Bad Gaia DR3 spectrophotometry! Not saving."
      endif else print, 'No Gaia DR3 spectrophotometry.'
   endif
endif

if keyword_set(gaiaonly) then begin
   if n_elements(sedfile) then free_lun, sedlun
   return
endif

; get GALEX
;print, 'Querying Bianchi+ (2011) for GALEX DR5 FUV and NUV...'
;qgalex=QueryVizier('II/312/ais',star,dist/60.,/silent,/all,cfa=cfa)
;if n_elements(qgalex) gt 1 then begin ; multiple matches. Just take the closest.
;  index = where(qgalex[*]._r eq min(qgalex[*]._r))
;  qgalex = qgalex[index]
;endif
;if long(tag_exist(qgalex,'fuv_6',/quiet)) ne 0L then begin
;  printf, sedlun, '# GALEX DR5 (Bianchi+2011; II/312/ais)' 
;  if qgalex.fuv_6 gt -99 and not finite(qgalex.e_fuv_6,/nan) then printf,sedlun,'GALEX_GALEX.FUV',qgalex.fuv_6+18.82-0.07,max([0.1,qgalex.e_fuv_6],/nan),qgalex.e_fuv_6 
;  if qgalex.nuv_6 gt -99 and not finite(qgalex.e_nuv_6,/nan) then printf,sedlun,'GALEX_GALEX.NUV',qgalex.nuv_6+20.08-0.07,max([0.1,qgalex.e_nuv_6],/nan),qgalex.e_nuv_6
;endif else print, "No GALEX match."
print, 'Querying Bianchi+ (2017; II/335/galex_ais) for GALEX GR7 FUV and NUV...'
qgalex=QueryVizier('II/335/galex_ais',star,dist/60.,/silent,/all,cfa=cfa)
if n_elements(qgalex) gt 1 then begin ; multiple matches. Just take the closest.
  index = where(qgalex[*]._r eq min(qgalex[*]._r))
  qgalex = qgalex[index]
endif
if long(tag_exist(qgalex,'fuv_6',/quiet)) ne 0L then begin
  printf, sedlun, '# GALEX GR6+7 (Bianchi+2017; II/335/galex_ais)' 
  if qgalex.fuv gt -99 and not finite(qgalex.e_fuv,/nan) then printf,sedlun,'GALEX_GALEX.FUV',qgalex.fuv,max([0.1,qgalex.e_fuv],/nan),qgalex.e_fuv
  if qgalex.nuv gt -99 and not finite(qgalex.e_nuv,/nan) then printf,sedlun,'GALEX_GALEX.NUV',qgalex.nuv,max([0.1,qgalex.e_nuv],/nan),qgalex.e_nuv
endif else print, "No GALEX match."

; get Mermilliod 1991 UBV
print, 'Querying Mermilliod 1991 UBV catalog...'
qmermilliod91=QueryVizier('II/168/ubvmeans',star,dist/60.,/silent,/all,cfa=cfa)
if n_elements(qmermilliod91) gt 1 then begin ; multiple entries. Find the one that explicitly matches input name
  index = where(strcmp(qmermilliod91(*).simbadname,star))
  qmermilliod91 = qmermilliod91(index)
endif
if long(tag_exist(qmermilliod91,'vmag',/quiet)) ne 0L then begin
  if qmermilliod91.u_b+qmermilliod91.b_v+qmermilliod91.vmag gt -9 and not finite(qmermilliod91.e_u_b,/nan) and not finite(qmermilliod91.e_b_v,/nan) $
  and not finite(qmermilliod91.e_vmag,/nan) then begin
    printf, sedlun, '# Mermilliod (1991); II/168/ubvmeans'
    printf, sedlun, 'OHP_Cam120.U',qmermilliod91.u_b+qmermilliod91.b_v+qmermilliod91.vmag,max([0.01,sqrt(qmermilliod91.e_u_b^2+qmermilliod91.e_b_v^2+qmermilliod91.e_vmag^2)],/nan)
  endif
  if qmermilliod91.b_v+qmermilliod91.vmag gt -9 and not finite(qmermilliod91.e_b_v,/nan) and not finite(qmermilliod91.e_vmag,/nan) then printf,sedlun,'OHP_Cam120.B',qmermilliod91.b_v+qmermilliod91.vmag,max([0.01,sqrt(qmermilliod91.e_b_v^2+qmermilliod91.e_vmag^2)],/nan)
  if qmermilliod91.vmag gt -9 and not finite(qmermilliod91.e_vmag,/nan) then printf,sedlun,'OHP_Cam120.V',qmermilliod91.vmag,max([0.01,qmermilliod91.e_vmag],/nan)
endif else print, "No Mermilliod 1991 match."

; get Tycho-2
print, 'Querying Hog+ (2000) for Tycho-2 B_T and VT...'
qtyc2=QueryVizier('I/259/tyc2',star,dist/60.,/silent,/all,cfa=cfa)
;if n_elements(qtyc2) gt 1 then begin ; XXXX:YYYY:1-2 (blend?), so take -1
;  index = where(qtyc2[*].tyc3 eq 1)
;  qtyc2 = qtyc2[index]
;endif
if n_elements(qtyc2) gt 1 then begin ; multiple matches, so take closest
  index = where(qtyc2[*]._r eq min(qtyc2[*]._r))
  qtyc2 = qtyc2[index]
endif
if long(tag_exist(qtyc2,'BTmag',/quiet)) ne 0L then begin
  printf, sedlun, '# Tycho-2 (Hog+ 2000; I/259/tyc2)'
  if qtyc2.btmag gt -9 and not finite(qtyc2.e_btmag,/nan) then printf,sedlun,'TYCHO_TYCHO.B_MvB',qtyc2.btmag,max([0.01,qtyc2.e_btmag],/nan),qtyc2.e_btmag
  if qtyc2.vtmag gt -9 and not finite(qtyc2.e_vtmag,/nan) then printf,sedlun,'TYCHO_TYCHO.V_MvB',qtyc2.vtmag,max([0.01,qtyc2.e_vtmag],/nan),qtyc2.e_vtmag
endif else begin
	print, "No Tycho-2 match."
	qtyc2={btmag:-99.,vtmag:-99.}
endelse

; get 2MASS
print, 'Querying Cutri+ (2003) 2MASS catalog for JHK...'
qtmass=QueryVizier('II/246/out',star,dist/60.,/silent,/all,cfa=cfa)
if n_elements(qtmass) gt 1 then begin
  index = where(qtmass[*]._r eq min(qtmass[*]._r))
  qtmass = qtmass[index]
endif
if long(tag_exist(qtmass,'Hmag',/quiet)) ne 0L then begin
  printf, sedlun, '# 2MASS (Cutri+2003; II/246/out)'
  if qtmass.Jmag gt -9 then printf,sedlun,'2MASS_2MASS.J',qtmass.Jmag,max([0.02,qtmass.e_Jmag],/nan),qtmass.e_Jmag
  if qtmass.Hmag gt -9 then printf,sedlun,'2MASS_2MASS.H',qtmass.Hmag,max([0.02,qtmass.e_Hmag],/nan),qtmass.e_Hmag ; changed max from 0.01->0.02 like mkticsed.pro,
  if qtmass.Kmag gt -9 then printf,sedlun,'2MASS_2MASS.Ks',qtmass.Kmag,max([0.02,qtmass.e_Kmag],/nan),qtmass.e_Kmag ; not sure why
endif else print, "No 2MASS match."


; get Paunzen 2015 Stromgren
print, 'Querying Paunzen+ (2015) for Stromgren-Crawford ubvy...'
qpaunzen15=QueryVizier('J/A+A/580/A23/catalog',star,dist/60.,/silent,/all,cfa=cfa)
if long(tag_exist(qpaunzen15,'vmag',/quiet)) ne 0L then begin
  printf,sedlun,'# Stromgren ubvy (Paunzen+2015; J/A+A/580/A23/catalog)'
  ubvymags=str_conv(qpaunzen15.vmag,max([0.01,qpaunzen15.e_vmag]),qpaunzen15.b_y,max([0.01,qpaunzen15.e_b_y]),qpaunzen15.m1,max([0.01,qpaunzen15.e_m1]),qpaunzen15.c1,max([0.01,qpaunzen15.e_c1]))
  if ubvymags(0) gt -9 then printf,sedlun,'Generic_Stromgren.u',ubvymags(0),max([0.02,ubvymags(1)],/nan)
  if ubvymags(2) gt -9 then printf,sedlun,'Generic_Stromgren.v',ubvymags(2),max([0.02,ubvymags(3)],/nan); changed max from 0.01->0.02 like mkticsed.pro,
  if ubvymags(4) gt -9 then printf,sedlun,'Generic_Stromgren.b',ubvymags(4),max([0.02,ubvymags(5)],/nan); not sure why
  if ubvymags(6) gt -9 then printf,sedlun,'Generic_Stromgren.y',ubvymags(6),max([0.02,ubvymags(7)],/nan)
endif else print, 'No Paunzen+ (2015) match.'

; get ucac4 - check if APASS BV is actually just TYC2 BV, and adjust errors from integer to centimag
;print, 'Querying UCAC4/APASS (Zacharias+2012; I/322A/out) for optical...'
;qucac4=QueryVizier('UCAC4',star,dist/60.,/silent,/all,cfa=cfa)
;if n_elements(qucac4) gt 1 then begin ; multiple matches. Just take the closest.
;  index = where(qucac4[*]._r eq min(qucac4[*]._r))
;  qucac4 = qucac4[index]
;endif
;if long(tag_exist(qucac4,'bmag',/quiet)) ne 0L then begin
;  printf,sedlun,'; UCAC4 (Zacharias+2012; I/322A/out)'
;  if qucac4.bmag ne qtyc2.btmag and qucac4.bmag gt -9 and qucac4.e_bmag ne 99 then printf,sedlun,'; B',qucac4.bmag,max([0.01,qucac4.e_bmag*0.01],/nan),qucac4.e_bmag*0.01
;  if qucac4.vmag ne qtyc2.vtmag and qucac4.bmag gt -9 and qucac4.e_vmag ne 99 then printf,sedlun,'; V',qucac4.vmag,max([0.01,qucac4.e_vmag*0.01],/nan),qucac4.e_vmag*0.01
;  if qucac4.gmag gt -9 then printf,sedlun,'; gSDSS',qucac4.gmag,max([0.01,qucac4.e_gmag*0.01],/nan),qucac4.e_gmag*0.01
;  if qucac4.rmag gt -9 then printf,sedlun,'; rSDSS',qucac4.rmag,max([0.01,qucac4.e_rmag*0.01],/nan),qucac4.e_rmag*0.01
;  if qucac4.imag gt -9 then printf,sedlun,'; iSDSS',qucac4.imag,max([0.01,qucac4.e_imag*0.01],/nan),qucac4.e_imag*0.01
;  if qucac4.Jmag gt -9 then printf,sedlun,'; J2M',qucac4.Jmag,max([0.01,qucac4.e_Jmag],/nan),qucac4.e_Jmag,' ; from UCAC4 catalog'
;  if qucac4.Hmag gt -9 then printf,sedlun,'; H2M',qucac4.Hmag,max([0.01,qucac4.e_Hmag],/nan),qucac4.e_Hmag
;  if qucac4.Kmag gt -9 then printf,sedlun,'; K2M',qucac4.Kmag,max([0.01,qucac4.e_Kmag],/nan),qucac4.e_Kmag
;endif else print, "No UCAC4 match."

;qapass=QueryVizier('APASS',star,dist/60.,/silent,/all,cfa=cfa)
;if n_elements(qapass) gt 1 then begin ; multiple matches. Just take the closest.
;  index = where(qapass[*]._r eq min(qapass[*]._r))
;  qapass = qapass[index]
;endif
;if long(tag_exist(qapass,'vmag',/quiet)) ne 0L then begin
;  bvgrimag=str_conv(qapass.bmag,max([0.01,qapass.e_bmag]),qapass.vmag,max([0.01,qapass.e_vmag]),qapass.g\'mag,max([0.01,qapass.e_g\'mag]),qapass.r\'mag,max([0.01,qapass.e_r\'mag]),qapass.i\'mag,max([0.01,qapass.e_i\'mag]))
;  if bvgrimag(0) gt -9 then printf,sedlun,'B',bvgrimag(0),max([0.01,bvgrimag(1)],/nan),qapass.e_bmag
;  if bvgrimag(2) gt -9 then printf,sedlun,'V',bvgrimag(2),max([0.01,bvgrimag(4)],/nan),qapass.e_vmag
;  if bvgrimag(4) gt -9 then printf,sedlun,'gSDSS',bvgrimag(4),max([0.01,bvgrimag(5)],/nan),qapass.e_g\'mag
;  if bvgrimag(6) gt -9 then printf,sedlun,'rSDSS',bvgrimag(6),max([0.01,bvgrimag(7)],/nan),qapass.e_r\'mag
;  if bvgrimag(8) gt -9 then printf,sedlun,'iSDSS',bvgrimag(8),max([0.01,bvgrimag(9)],/nan),qapass.e_i\'mag
;endif else print, "No APASS match."


; get ALLWISE
print, 'Querying AllWISE (Cutri+2013; II/328/allwise) for IR WISE1-4...'
qwise=QueryVizier('II/328/allwise',star,dist/60.,/silent,/all,cfa=cfa)
if n_elements(qwise) gt 1 then begin ; multiple matches. Just take the closest.
  index = where(qwise[*]._r eq min(qwise[*]._r))
  qwise = qwise[index]
endif
if long(tag_exist(qwise,'w1mag',/quiet)) ne 0L then begin
  printf,sedlun,'# AllWISE (Cutri+2013; II/328/allwise)' 
  if qwise.w1mag gt -9 and not finite(qwise.e_w1mag,/nan) then printf,sedlun,'WISE_WISE.W1',qwise.w1mag,max([0.03d,qwise.e_w1mag],/nan),qwise.e_w1mag
  if qwise.w2mag gt -9 and not finite(qwise.e_w2mag,/nan) then printf,sedlun,'WISE_WISE.W2',qwise.w2mag,max([0.03d,qwise.e_w2mag],/nan),qwise.e_w2mag ; changed 0.01->0.03d like EXOFASTv2 mkticsed.pro,
  if qwise.w3mag gt -9 and not finite(qwise.e_w3mag,/nan) then printf,sedlun,'WISE_WISE.W3',qwise.w3mag,max([0.03d,qwise.e_w3mag],/nan),qwise.e_w3mag ; but not clear why.
  if qwise.w4mag gt -9 and not finite(qwise.e_w4mag,/nan) then printf,sedlun,'WISE_WISE.W4',qwise.w4mag,max([0.10d,qwise.e_w4mag],/nan),qwise.e_w4mag
endif else print, "No AllWISE match."

print, "Querying Abdurro'uf+ (2022) for APOGEE-2 DR17 RVs..."
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
endif else print, "No APOGEE-2 DR17 RV match."


if n_elements(sedfile) then free_lun, sedlun
end
