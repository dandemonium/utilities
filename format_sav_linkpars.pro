pro format_sav_linkpars, mcmcfile, savname, mist=mist, circ=circ, ellip=ellip, beam=beam, reflect=reflect
restore, mcmcfile
chi2 = *(mcmcss.chi2)
chi2 = reform(chi2,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
burnndx = getburnndx(chi2,goodchains=goodchains)
ntran = mcmcss.ntran
nband = mcmcss.nband
nstars = mcmcss.nstars

sedfile = mcmcss.sedfile
sedrange = mcmcss.sedrange
gravitysun = mcmcss.constants.gravitysun
teffemfloor = mcmcss.teffemfloor
fehfloor = mcmcss.fehemfloor
rstarfloor = mcmcss.rstaremfloor
agefloor = mcmcss.ageemfloor
emrange = mcmcss.emrange

if keyword_set(mist) then begin
	par1d = mcmcss.star[*].age
	par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
	age = par_2d[*, burnndx:-1,goodchains]
	age = reform(age,npars[1]*npars[2], nstars)

	par1d = mcmcss.star[*].eep
	par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
	eep = par_2d[*, burnndx:-1,goodchains]
	eep = reform(eep, npars[1]*npars[2], nstars)

 	par1d = mcmcss.star[*].initfeh
	par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
	initfeh = par_2d[*, burnndx:-1,goodchains]
	initfeh = reform(initfeh, npars[1]*npars[2], nstars)
endif
par1d = mcmcss.band[*].thermal
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
thermal = par_2d[burnndx:-1,goodchains]
npars = size(thermal)
thermal = reform(thermal,npars[1]*npars[2]);,mcmcss.nsteps*mcmcss.nchains)

;for i=0, nstars-1 do begin
par1d = mcmcss.star[*].rstar
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
rstar = par_2d[*, burnndx:-1,goodchains]
rstar = reform(rstar, npars[1]*npars[2], nstars)
par1d = mcmcss.star[*].rstarsed
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
rstarsed = par_2d[*, burnndx:-1,goodchains]
rstarsed = reform(rstarsed, npars[1]*npars[2], nstars)

par1d = mcmcss.star[*].lstar
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
lstar = par_2d[*, burnndx:-1,goodchains]
lstar = reform(lstar, npars[1]*npars[2], nstars)

par1d = mcmcss.star[*].logmstar
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
logmstar = par_2d[*, burnndx:-1,goodchains]
logmstar = reform(logmstar, npars[1]*npars[2], nstars)

par1d = mcmcss.star[*].logg 
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
logg = par_2d[*, burnndx:-1,goodchains]
logg = reform(logg, npars[1]*npars[2], nstars)

par1d = mcmcss.star[*].feh 
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
feh = par_2d[*, burnndx:-1,goodchains]
feh = reform(feh, npars[1]*npars[2], nstars)

par1d = mcmcss.star[*].teff
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
teff = par_2d[*, burnndx:-1,goodchains]
teff = reform(teff, npars[1]*npars[2], nstars)
par1d = mcmcss.star[*].teffsed
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
teffsed = par_2d[*, burnndx:-1,goodchains]
teffsed = reform(teffsed, npars[1]*npars[2], nstars)

par1d = mcmcss.star[*].av 
par_2d = (reform(par1d.value, nstars, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
av = par_2d[*, burnndx:-1,goodchains]
av = reform(av,npars[1]*npars[2], nstars)

lstarsed = 4d0*!dpi*rstarsed^2*teffsed^4*mcmcss.constants.sigmab/mcmcss.constants.lsun*mcmcss.constants.rsun^2 ;; lsun

par1d = mcmcss.star[0].errscale
par_2d = (reform(par1d.value, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
errscale = par_2d[burnndx:-1,goodchains]
errscale = reform(errscale, npars[1]*npars[2])

;endfor
;par1d = mcmcss.star[0].logmstar
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;logmstar = par_2d[burnndx:-1,goodchains]
;logmstar = reform(logmstar,npars[1]*npars[2])

;par1d = mcmcss.star[0].lstar
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;lstar1 = par_2d[burnndx:-1,goodchains]
;lstar1 = reform(lstar1,npars[1]*npars[2])

;par1d = mcmcss.star[1].lstar
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;lstar2 = par_2d[burnndx:-1,goodchains]
;star2 = reform(lstar2,npars[1]*npars[2])

;par1d = mcmcss.star[1].logmstar
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;logmstar2 = par_2d[burnndx:-1,goodchains]
;logmstar2 = reform(logmstar2,npars[1]*npars[2])

;par1d = mcmcss.star[0].mstar
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;mstar = par_2d[burnndx:-1,goodchains]
;mstar = reform(mstar,npars[1]*npars[2])

;par1d = mcmcss.star[0].logg 
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;loggstar = par_2d[burnndx:-1,goodchains]
;loggstar = reform(loggstar,npars[1]*npars[2])

;par1d = mcmcss.star[0].feh 
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;feh = par_2d[burnndx:-1,goodchains]
;feh = reform(feh,npars[1]*npars[2])

;par1d = mcmcss.star[1].feh
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;fehp = par_2d[burnndx:-1,goodchains]
;fehp = reform(fehp,npars[1]*npars[2])


;par1d = mcmcss.star[0].teff 
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;teff = par_2d[burnndx:-1,goodchains]
;teff = reform(teff,npars[1]*npars[2])

;par1d = mcmcss.star[1].teff
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;teffp = par_2d[burnndx:-1,goodchains]
;teffp = reform(teffp,npars[1]*npars[2])

;par1d = mcmcss.planet.loggp 
;par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
;loggp = par_2d[burnndx:-1,goodchains]
;loggp = reform(loggp,npars[1]*npars[2])



par1d = mcmcss.planet.p
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
rprs = par_2d[burnndx:-1,goodchains]
rprs = reform(rprs,npars[1]*npars[2])

par1d = mcmcss.planet.rpsun
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
rp = par_2d[burnndx:-1,goodchains]
rp = reform(rp,npars[1]*npars[2])

par1d = mcmcss.planet.mpsun
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
mp = par_2d[burnndx:-1,goodchains]
mp = reform(mp,npars[1]*npars[2])

par1d = mcmcss.planet.ar
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
ar = par_2d[burnndx:-1,goodchains]
ar = reform(ar, npars[1]*npars[2])
par1d = mcmcss.planet.a
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
a = par_2d[burnndx:-1,goodchains]
a = reform(a, npars[1]*npars[2])

par1d = mcmcss.planet.esinw
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
esinw = par_2d[burnndx:-1,goodchains]
esinw= reform(esinw,npars[1]*npars[2])

par1d = mcmcss.planet.ecosw
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
ecosw = par_2d[burnndx:-1,goodchains]
ecosw= reform(ecosw,npars[1]*npars[2])

par1d = mcmcss.planet.e
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
e = par_2d[burnndx:-1,goodchains]
e= reform(e,npars[1]*npars[2])
par1d = mcmcss.planet.omegadeg
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
omegadeg = par_2d[burnndx:-1,goodchains]
omegadeg= reform(omegadeg,npars[1]*npars[2])


par1d = mcmcss.planet.cosi
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
cosi = par_2d[burnndx:-1,goodchains]
cosi= reform(cosi,npars[1]*npars[2])
par1d = mcmcss.planet.ideg
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
ideg = par_2d[burnndx:-1,goodchains]
ideg= reform(ideg,npars[1]*npars[2])

par1d = mcmcss.planet.tc
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
tc = par_2d[burnndx:-1,goodchains]
tc= reform(tc,npars[1]*npars[2])

par1d = mcmcss.planet.ts
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
ts = par_2d[burnndx:-1,goodchains]
ts = reform(ts,npars[1]*npars[2])
par1d = mcmcss.planet.tp
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
tp = par_2d[burnndx:-1,goodchains]
tp = reform(tp,npars[1]*npars[2])

par1d = mcmcss.planet.period
par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
period = par_2d[burnndx:-1,goodchains]
period= reform(period,npars[1]*npars[2])

if keyword_set(ellip) then begin
	par1d = mcmcss.band[*].ellip
	par_2d = (reform(par1d.value, nband, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
	ellip = par_2d[*, burnndx:-1,goodchains]
	ellip = reform(ellip, npars[1]*npars[2], nband)
endif
if keyword_set(reflect) then begin
	par1d = mcmcss.band[*].reflect
	par_2d = (reform(par1d.value, nband, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
	reflect = par_2d[*, burnndx:-1,goodchains]
	reflect = reform(reflect, npars[1]*npars[2], nband)
endif
if keyword_set(beam) then begin
	par1d = mcmcss.planet.beam
	par_2d = (reform(par1d.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
	beam = par_2d[burnndx:-1,goodchains]
	beam = reform(beam,npars[1]*npars[2])
endif

if keyword_set(ntran) then begin
	par1d = mcmcss.transit[*].f0
	par_2d = (reform(par1d.value, ntran, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
	f0 = par_2d[*, burnndx:-1,goodchains]
	f0 = reform(f0,npars[1]*npars[2], ntran)
	par1d = mcmcss.transit[*].dilute
	par_2d = (reform(par1d.value, ntran, mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))
	dilute = par_2d[*, burnndx:-1,goodchains]
	dilute = reform(dilute,npars[1]*npars[2], ntran)
 endif

;burnndx = mcmcss.burnndx;
undefine, mcmcss, chi2, burnndx, par1d, par_2d
save,/variables,filename = savname 
;writecol, savname+'.print', thermal, rstar, loggstar, feh, fehp, av, teff, teffp, loggp, rprs, rp
;forprint, thermal, rstar, loggstar, feh, fehp, av, teff, teffp, loggp, rprs, rp, /nocomment, textout=savname+'.print',format='(2(f0.6,x))'
return
end
