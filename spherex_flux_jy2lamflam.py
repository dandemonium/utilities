"""Convert SPHEREx spectrophotometry from IRSA:
	lam(um)	lam_width(um)	Flux(uJy)	Flux_err(uJy)
	into format expected by EXOFASTv2:
	lam(um)	Flux(cgs)	Flux_err(cgs)
"""

import numpy as np, sys

def convert_to_lamflam(lam, f_ujy):
	""" Converts input flux F in micro-Janksy
		to F_lambda, then to lambda*F_lambda,
		where lambda is in microns.
		Formula from Appendix 23 of kaSTScI's
		NICMOS Instrument Handbook for Cycle
		11, version 5:
		flam = beta * fnu / lam"""

	beta = 3e-9 #to get erg/s/cm^2/um
	return np.float64((beta*f_ujy*1e-6)/lam)
	
fname = sys.argv[1]
star = sys.argv[2]
lam, width, flux_ujy, flux_err_ujy = np.loadtxt(fname, unpack=True, comments='#')
f_cgs = convert_to_lamflam(lam, flux_ujy)
f_err_cgs = convert_to_lamflam(lam, flux_err_ujy)
np.savetxt(star+".SPHEREx.sed", np.array([lam,f_cgs,f_err_cgs]).transpose(), fmt=('%.8f','%.8e','%.8e'), delimiter='\t')
