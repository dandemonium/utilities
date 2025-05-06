import numpy as np, sys, argparse as ap

def eggleton(q):
	"""Estimate scaled Roche radius for given
	   orbital separation a (in Rsun) and
	   mass ratio q, per Eggleton (1983, ApJ)."""
	   
	qthird = q**(-1./3.)
	num = 0.49 * qthird * qthird
	denom = 0.6 * qthird * qthird
	denom += np.log(1. + qthird)
	   
	return denom / num # a/R_roche
	
def paczynski(q):
	"""Estimate scaled Roche radius for given
	   orbital separation a (in Rsun) and
	   mass ratio q, per Paczynski (1971, ARAA)."""
	if (q >= 20) or (q <= 0):
		print("Mass ratio is out of bounds (0 < q < 20).")
		return
	f1 = 0.38 + 0.2*np.log(1./q)
	f2 = 0.46224 * (1. + q)**(-1./3.)
	return 1. / max(f1,f2) # a/R_roche

def kepler(q, p, m1):
	"""Use Kepler's 3rd Law to compute semimajor axis."""
	G = 6.67e-11
	msun = 1.9891e30
	rsun = 6.957e8
	return ((p * 24. * 3600.)**2. * G * m1 * (1. + q) * msun / (4 * np.pi**2.))**(1./3.) / rsun # a, in Rsun
	
parser = ap.ArgumentParser(prog='roche', description="Estimate Roche radius for given separation, mass ratio, and method.")
parser.add_argument("-q", "--q", help="Mass ratio, q (0 < q < 1)", type=float, default=1.0)
parser.add_argument("-a", "--a", help="Orbital separation, a (Rsun)", type=float)
parser.add_argument("-p", "--period", help="Orbital period, P (days)", type=float)
parser.add_argument("-r", "--radius", help="Stellar radius, R (Rsun)", type=float, default=1.0)
parser.add_argument("-m", "--mass", help="Primary stellar mass, M (Msun))", type=float, default=1.0)
args = parser.parse_args()

if args.q == None:
	q = 1.0
else: q = args.q
roche_eggleton = eggleton(q)
roche_paczynski = paczynski(q)

print("Eggleton (1983, ApJ) a/R: ", roche_eggleton)
print("Paczynski (1971, ARAA) a/R: ", roche_paczynski)

# compare scaled Roche radius to scaled stellar radius
r = args.radius
print("Setting R_1 = %f R_sun for the following calculations:" % r)
if args.a != None:
	print("User-provided scaled orbital separation a/R* is ", args.a/r)
elif args.period != None:
	m = args.mass
	p = args.period
	print("Scaled orbital separation a/R* from Kepler's 3rd Law is ", kepler(q, p, m)/r)

