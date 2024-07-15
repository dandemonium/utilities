import numpy as np, sys

filename = sys.argv[1]

teff, eteff, logg, elogg, feh, efeh = np.loadtxt(filename, unpack=True)

eteff_wt = np.sum(eteff**-2.)
teff_wt = np.sum(teff/(eteff**2.)) / eteff_wt
elogg_wt = np.sum(elogg**-2.) 
logg_wt = np.sum(logg/(elogg**2.)) / elogg_wt
efeh_wt = np.sum(efeh**-2.)
feh_wt = np.sum(feh/(efeh**2.)) / efeh_wt

print("Error-weighted Teff, log(g), and [Fe/H] are:")
print("teff_0 {} {}".format(teff_wt, 1.0/np.sqrt(eteff_wt)))
print("logg_0 {} # {}".format(logg_wt, 1.0/np.sqrt(elogg_wt)))
print("feh_0 {} {}".format(feh_wt, 1.0/np.sqrt(efeh_wt)))
print("feh_1 feh_0 0.15")
