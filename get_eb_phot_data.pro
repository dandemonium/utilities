pro	get_eb_phot_data,star,outfile,gaiaxp=gaiaxp
lun =1
if n_params() lt 2 then begin
  print,'syntax: get_eb_phot_data,star,outfile,gaiaxp=gaiaxp'
  retall
endif
close,lun
openw,lun,outfile

; get GALEX
print, 'Querying Bianchi+ (2011) for GALEX DR5 FUV and NUV...'
qgalex=QueryVizier('II/312/ais',star,6./60.,/silent,/all,/CFA)
if n_elements(qgalex) gt 1 then begin ; multiple matches. Just take the closest.
  index = where(qgalex[*]._r eq min(qgalex[*]._r))
  qgalex = qgalex[index]
endif
if long(tag_exist(qgalex,'fuv_6',/quiet)) ne 0L then begin
  if qgalex.fuv_6 gt -99 and not finite(qgalex.e_fuv_6,/nan) then printf,lun,'galFUV',qgalex.fuv_6,max([0.1,qgalex.e_fuv_6],/nan),qgalex.e_fuv_6, '# GALEX DR5 (Bianchi+2011; II/312/ais)'  
  if qgalex.nuv_6 gt -99 and not finite(qgalex.e_nuv_6,/nan) then printf,lun,'galNUV',qgalex.nuv_6,max([0.1,qgalex.e_nuv_6],/nan),qgalex.e_nuv_6
endif else print, "No GALEX match."

; get Mermilliod 1991 UBV
print, 'Querying Mermilliod 1991 UBV catalog...'
qmermilliod91=QueryVizier('II/168/ubvmeans',star,6./60.,/silent,/all,/CFA)
if n_elements(qmermilliod91) gt 1 then begin ; multiple entries. Find the one that explicitly matches input name
  index = where(strcmp(qmermilliod91(*).simbadname,star))
  qmermilliod91 = qmermilliod91(index)
endif
if long(tag_exist(qmermilliod91,'vmag',/quiet)) ne 0L then begin
  if qmermilliod91.u_b+qmermilliod91.b_v+qmermilliod91.vmag gt -9 and not finite(qmermilliod91.e_u_b,/nan) and not finite(qmermilliod91.e_b_v,/nan) and not finite(qmermilliod91.e_vmag,/nan) then printf,lun,'U',qmermilliod91.u_b+qmermilliod91.b_v+qmermilliod91.vmag,max([0.01,sqrt(qmermilliod91.e_u_b^2+qmermilliod91.e_b_v^2+qmermilliod91.e_vmag^2)],/nan),'# Mermilliod (1991); II/168/ubvmeans'
  if qmermilliod91.b_v+qmermilliod91.vmag gt -9 and not finite(qmermilliod91.e_b_v,/nan) and not finite(qmermilliod91.e_vmag,/nan) then printf,lun,'B',qmermilliod91.b_v+qmermilliod91.vmag,max([0.01,sqrt(qmermilliod91.e_b_v^2+qmermilliod91.e_vmag^2)],/nan)
  if qmermilliod91.vmag gt -9 and not finite(qmermilliod91.e_vmag,/nan) then printf,lun,'V',qmermilliod91.vmag,max([0.01,qmermilliod91.e_vmag],/nan)
endif else print, "No Mermilliod 1991 match."

; get Tycho-2
print, 'Querying Hog+ (2000) for Tycho-2 B_T and VT...'
qtyc2=QueryVizier('I/259/tyc2',star,6./60.,/silent,/all,/CFA)
;if n_elements(qtyc2) gt 1 then begin ; XXXX:YYYY:1-2 (blend?), so take -1
;  index = where(qtyc2[*].tyc3 eq 1)
;  qtyc2 = qtyc2[index]
;endif
if n_elements(qtyc2) gt 1 then begin ; multiple matches, so take closest
  index = where(qtyc2[*]._r eq min(qtyc2[*]._r))
  qtyc2 = qtyc2[index]
endif
if long(tag_exist(qtyc2,'BTmag',/quiet)) ne 0L then begin
  if qtyc2.btmag gt -9 and not finite(qtyc2.e_btmag,/nan) then printf,lun,'BT',qtyc2.btmag,max([0.01,qtyc2.e_btmag],/nan),qtyc2.e_btmag,'# Tycho-2 (Hog+ 2000; I/259/tyc2)'
  if qtyc2.vtmag gt -9 and not finite(qtyc2.e_vtmag,/nan) then printf,lun,'VT',qtyc2.vtmag,max([0.01,qtyc2.e_vtmag],/nan),qtyc2.e_vtmag
endif else begin
	print, "No Tycho-2 match."
	qtyc2={btmag:-99.,vtmag:-99.}
endelse

; get 2MASS
print, 'Querying Cutri+ (2003) 2MASS catalog for JHK...'
qtmass=QueryVizier('II/246/out',star,6./60.,/silent,/all,/CFA)
if n_elements(qtmass) gt 1 then begin
  index = where(qtmass[*]._r eq min(qtmass[*]._r))
  qtmass = qtmass[index]
endif
if long(tag_exist(qtmass,'Hmag',/quiet)) ne 0L then begin
  if qtmass.Jmag gt -9 then printf,lun,'J2M',qtmass.Jmag,max([0.02,qtmass.e_Jmag],/nan),qtmass.e_Jmag, '# 2MASS (Cutri+2003; II/246/out)'
  if qtmass.Hmag gt -9 then printf,lun,'H2M',qtmass.Hmag,max([0.02,qtmass.e_Hmag],/nan),qtmass.e_Hmag ; changed max from 0.01->0.02 like mkticsed.pro,
  if qtmass.Kmag gt -9 then printf,lun,'K2M',qtmass.Kmag,max([0.02,qtmass.e_Kmag],/nan),qtmass.e_Kmag ; not sure why
endif else print, "No 2MASS match."


; get Paunzen 2015 Stromgren
print, 'Querying Paunzen+ (2015) for Stromgren-Crawford ubvy...'
qpaunzen15=QueryVizier('J/A+A/580/A23/catalog',star,6./60.,/silent,/all,/CFA)
if long(tag_exist(qpaunzen15,'vmag',/quiet)) ne 0L then begin
  ubvymags=str_conv(qpaunzen15.vmag,max([0.01,qpaunzen15.e_vmag]),qpaunzen15.b_y,max([0.01,qpaunzen15.e_b_y]),qpaunzen15.m1,max([0.01,qpaunzen15.e_m1]),qpaunzen15.c1,max([0.01,qpaunzen15.e_c1]))
  if ubvymags(0) gt -9 then printf,lun,'uStr',ubvymags(0),max([0.02,ubvymags(1)],/nan),'# Stromgren ubvy (Paunzen+2015; J/A+A/580/A23/catalog)'
  if ubvymags(2) gt -9 then printf,lun,'vStr',ubvymags(2),max([0.02,ubvymags(3)],/nan); changed max from 0.01->0.02 like mkticsed.pro,
  if ubvymags(4) gt -9 then printf,lun,'bStr',ubvymags(4),max([0.02,ubvymags(5)],/nan); not sure why
  if ubvymags(6) gt -9 then printf,lun,'yStr',ubvymags(6),max([0.02,ubvymags(7)],/nan)
endif else print, 'No Paunzen+ (2015) match.'

; get ucac4 - check if APASS BV is actually just TYC2 BV, and adjust errors from integer to centimag
;print, 'Querying UCAC4/APASS (Zacharias+2012; I/322A/out) for optical...'
;qucac4=QueryVizier('UCAC4',star,6./60.,/silent,/all,/CFA)
;if n_elements(qucac4) gt 1 then begin ; multiple matches. Just take the closest.
;  index = where(qucac4[*]._r eq min(qucac4[*]._r))
;  qucac4 = qucac4[index]
;endif
;if long(tag_exist(qucac4,'bmag',/quiet)) ne 0L then begin
;  printf,lun,'; UCAC4 (Zacharias+2012; I/322A/out)'
;  if qucac4.bmag ne qtyc2.btmag and qucac4.bmag gt -9 and qucac4.e_bmag ne 99 then printf,lun,'; B',qucac4.bmag,max([0.01,qucac4.e_bmag*0.01],/nan),qucac4.e_bmag*0.01
;  if qucac4.vmag ne qtyc2.vtmag and qucac4.bmag gt -9 and qucac4.e_vmag ne 99 then printf,lun,'; V',qucac4.vmag,max([0.01,qucac4.e_vmag*0.01],/nan),qucac4.e_vmag*0.01
;  if qucac4.gmag gt -9 then printf,lun,'; gSDSS',qucac4.gmag,max([0.01,qucac4.e_gmag*0.01],/nan),qucac4.e_gmag*0.01
;  if qucac4.rmag gt -9 then printf,lun,'; rSDSS',qucac4.rmag,max([0.01,qucac4.e_rmag*0.01],/nan),qucac4.e_rmag*0.01
;  if qucac4.imag gt -9 then printf,lun,'; iSDSS',qucac4.imag,max([0.01,qucac4.e_imag*0.01],/nan),qucac4.e_imag*0.01
;  if qucac4.Jmag gt -9 then printf,lun,'; J2M',qucac4.Jmag,max([0.01,qucac4.e_Jmag],/nan),qucac4.e_Jmag,' ; from UCAC4 catalog'
;  if qucac4.Hmag gt -9 then printf,lun,'; H2M',qucac4.Hmag,max([0.01,qucac4.e_Hmag],/nan),qucac4.e_Hmag
;  if qucac4.Kmag gt -9 then printf,lun,'; K2M',qucac4.Kmag,max([0.01,qucac4.e_Kmag],/nan),qucac4.e_Kmag
;endif else print, "No UCAC4 match."

;qapass=QueryVizier('APASS',star,6./60.,/silent,/all,/CFA)
;if n_elements(qapass) gt 1 then begin ; multiple matches. Just take the closest.
;  index = where(qapass[*]._r eq min(qapass[*]._r))
;  qapass = qapass[index]
;endif
;if long(tag_exist(qapass,'vmag',/quiet)) ne 0L then begin
;  bvgrimag=str_conv(qapass.bmag,max([0.01,qapass.e_bmag]),qapass.vmag,max([0.01,qapass.e_vmag]),qapass.g\'mag,max([0.01,qapass.e_g\'mag]),qapass.r\'mag,max([0.01,qapass.e_r\'mag]),qapass.i\'mag,max([0.01,qapass.e_i\'mag]))
;  if bvgrimag(0) gt -9 then printf,lun,'B',bvgrimag(0),max([0.01,bvgrimag(1)],/nan),qapass.e_bmag
;  if bvgrimag(2) gt -9 then printf,lun,'V',bvgrimag(2),max([0.01,bvgrimag(4)],/nan),qapass.e_vmag
;  if bvgrimag(4) gt -9 then printf,lun,'gSDSS',bvgrimag(4),max([0.01,bvgrimag(5)],/nan),qapass.e_g\'mag
;  if bvgrimag(6) gt -9 then printf,lun,'rSDSS',bvgrimag(6),max([0.01,bvgrimag(7)],/nan),qapass.e_r\'mag
;  if bvgrimag(8) gt -9 then printf,lun,'iSDSS',bvgrimag(8),max([0.01,bvgrimag(9)],/nan),qapass.e_i\'mag
;endif else print, "No APASS match."


; get ALLWISE
print, 'Querying AllWISE (Cutri+2013; II/328/allwise) for IR WISE1-4...'
qwise=QueryVizier('II/328/allwise',star,6./60.,/silent,/all,/CFA)
if n_elements(qwise) gt 1 then begin ; multiple matches. Just take the closest.
  index = where(qwise[*]._r eq min(qwise[*]._r))
  qwise = qwise[index]
endif
if long(tag_exist(qwise,'w1mag',/quiet)) ne 0L then begin
  printf,lun,'# AllWISE (Cutri+2013; II/328/allwise)' 
  if qwise.w1mag gt -9 and not finite(qwise.e_w1mag,/nan) then printf,lun,'WISE1',qwise.w1mag,max([0.03d,qwise.e_w1mag],/nan),qwise.e_w1mag
  if qwise.w2mag gt -9 and not finite(qwise.e_w2mag,/nan) then printf,lun,'WISE2',qwise.w2mag,max([0.03d,qwise.e_w2mag],/nan),qwise.e_w2mag ; changed 0.01->0.03d like EXOFASTv2 mkticsed.pro,
  if qwise.w3mag gt -9 and not finite(qwise.e_w3mag,/nan) then printf,lun,'WISE3',qwise.w3mag,max([0.03d,qwise.e_w3mag],/nan),qwise.e_w3mag ; but not clear why.
  if qwise.w4mag gt -9 and not finite(qwise.e_w4mag,/nan) then printf,lun,'WISE4',qwise.w4mag,max([0.10d,qwise.e_w4mag],/nan),qwise.e_w4mag
endif else print, "No AllWISE match."

; get HIPPARCOS
;qhip=QueryVizier('I/311/hip2',star,10./60.,/silent,/all,/CFA)
;if long(tag_exist(qhip,'Plx',/quiet)) ne 0L then begin
;  if qhip.Plx gt -99 and not finite(qhip.e_Plx,/nan) then printf,lun,'; HIP_plx',qhip.Plx,qhip.e_Plx
;endif else print, "No Hipparcos match."

; get Gaia spectrophotometry
if keyword_set(gaiaxp) then begin
	gaiafile='gaia.bprp.'+STRJOIN(STRSPLIT(star, /EXTRACT), '')+'.dat'
	print, 'Querying Gaia DR3 Pt. 1 Main Src: BP/RP ext. cal. sampled mean spectrum (I/355/xpsample)'
	qgaia_bprp=QueryVizier('I/355/xpsample',star,3./60.,/silent,/all,/CFA)
	if n_elements(qgaia_bprp) gt 1 then begin ; multiple matches, so take closest
	  index = where(qgaia_bprp[*]._r eq min(qgaia_bprp[*]._r))
	  qgaia_bprp = qgaia_bprp[index]
	endif
	if long(tag_exist(qgaia_bprp,'flux',/quiet)) ne 0L then begin
		if qgaia_bprp[0].flux gt 0d0 and not finite(qgaia_bprp[0].e_flux,/nan) then begin	
            err_floor = dblarr(n_elements(qgaia_bprp.e_flux)) ; set 5% err floor, per K. Stassun (email)
			for i=0,n_elements(err_floor)-1 do begin
               err_floor[i] = max([qgaia_bprp[i].e_flux,0.05*qgaia_bprp[i].flux],/nan)
            endfor
			write_csv,gaiafile,qgaia_bprp.lambda,qgaia_bprp.flux,err_floor,qgaia_bprp.e_flux,table_header=[ $
			'#Gaia DR3 Pt. 1 Main Src: BP/RP ext. cal. sampled mean spectrum(2022 I/355/xpsample)'],$
			header=['#lambda(nm)','Flam(W/m^2/nm)','Flam_err_floor(5pct_W/m^2/nm)','Flam_err(W/m^2/nm)']
		endif else begin
		  	print, "Bad Gaia DR3 spectrophotometry! Not saving."
		endelse
	endif else print, 'No Gaia DR3 spectrophotometry.'
endif

close,lun
end
