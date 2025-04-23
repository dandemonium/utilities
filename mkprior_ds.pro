;+
; NAME:
;   STROM_CONV
; PURPOSE:
;    Translates stromgren color combinations from the catalog to
;    individual uvby magnitudes
; Modification 
;    2018-04-12: Jason Eastman, CfA
;                Renamed, documented, and cleaned up for distribution with EXOFASTv2
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
;+
; NAME:
;   mkprior_ds
; PURPOSE:
;    Grab key priors for EXOFASTv2 fits, including Av and Gaia zpt-corrected plx.
; Modification 
;    2018-04-12: Jason Eastman, CfA (as mkticsed.pro)
;                Renamed, documented, and cleaned up for distribution with EXOFASTv2
;    2025-04-23: Dan Stevens, UM-Duluth
;		  Switch from Gaia EDR3 to DR3; remove photometry queries;
;         add queries for APOGEE-2 DR17 spectroscopic parameters, Gaia DR3
;         RV parameters; rename mkprior_ds.pro
;-
pro mkprior_ds, ticid, priorfile=priorfile, france=france, ra=ra, dec=dec

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
      endif
   endfor
endif

if keyword_set(france) then cfa = 0B $
else cfa = 1B

;; match RA/Dec to TICv8 (closest)
if n_elements(ra) ne 0 and n_elements(dec) ne 0 then begin
   print, 'WARNING: querying by RA/Dec is less robust than querying by TIC ID and may lead to misidentification'   
   qtic = exofast_queryvizier('IV/38/tic',[ra,dec],2d0, /silent, /all, cfa=cfa)
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
   qtic = Exofast_Queryvizier('IV/38/tic','TIC ' + strtrim(ticid,2),/allcolumns,cfa=cfa)
endif else begin
   ;; query by supplied name (less robust)
   print, 'WARNING: querying by name is less robust than querying by TIC ID and may lead to misidentification'
   qtic = Exofast_Queryvizier('IV/38/tic',ticid,2d0,/allcolumns,cfa=cfa)
   qgaia = Exofast_Queryvizier('I/345/gaia2',ticid,2d0,/allcolumns,cfa=cfa)   

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
;ufeh = 0.08d0 > qtic.e__m_h_ 
;if finite(feh) and finite(ufeh) then begin
;   printf, priorlun, feh, '#',ufeh, format='("feh",x,f0.5,x,f0.5)'
;endif

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


;; use the Gaia ID to query the Gaia catalog
gaiaid = qtic.gaia
;; DR3 (print BP/RP/G, but leave it commented out)
qgaia3=Exofast_Queryvizier('I/355/gaiadr3',star,dist/60.,/silent,cfa=cfa,/all)
print, "Querying Gaia DR3..."
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
   endif
endif

print, "Querying Abdurro'uf+ (2022) for APOGEE-2 DR17 spectroscopic priors..."
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
endif else print, "No APOGEE-2 DR17 Teff, logg, [Fe/H] match."

print, "Querying Gaia DR3 (I/357/tbosb1, I/360/binmass) RV orbital solution for initial guesses (NOT YET SUPPORTED)..."
qgaiadr3_rv1 = QueryVizier('I/357/tbosb1', star, dist/60., /silent, cfa=cfa, /all)
qgaiadr3_rv2 = Queryvizier('I/360/binmass', star, dist/60., /silent, cfa=cfa, /all)

free_lun, priorlun

end

