#!/usr/bin/python

import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pacs import Pacsmap, _double_gauss, mirror
import argparse
from os.path import join
from christensen import datadir, psfdir, psf_vesta

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
parser.add_argument('-b', '--band', default="blue", choices=("blue", "red"),
                help="PACS band")
parser.add_argument('--binsize', default=0, type=float, help="bin size")
parser.add_argument('--rmax', default=0, type=int, help="radial profile extent")
args = parser.parse_args()

# calculate PSF radial profile
psf = Pacsmap(join(psfdir, psf_vesta[args.band]), comet=False,
                size=np.max((60, args.rmax)))
fov = psf.patch.shape[0]/2
psf.radprof(center=(fov,fov), rmax=args.rmax)
if args.debug:
    plt.imshow(psf.patch, origin="lower")
    plt.scatter(fov,fov)
    plt.show()
plt.scatter(psf.r, psf.rprof, marker='.', color=args.band)
ax = plt.gca()
ax.set_yscale('log')
if args.binsize:
    psf.radprof(binsize=args.binsize, center=(fov,fov), rmax=args.rmax)
    plt.errorbar(psf.r, psf.rprof, yerr=psf.rprof_e, fmt='x', color="g")
    for i in np.arange(0, args.rmax, args.binsize):
        plt.axvline(x=i, linestyle='--')
    plt.scatter(mirror(psf.r, sign=-1), mirror(psf.rprof), marker='x',
            color="g")
    np.savetxt(join(datadir, 'ascii',
        'PSF_{0}_{1}_prof.dat'.format(args.band, args.binsize)),
        np.transpose((psf.r, psf.rprof, psf.rprof_e)))
else:
    mask = [psf.rprof > 1e-6]
    np.savetxt(join(datadir, 'ascii',
        'PSF_{0}_{1}_prof.dat'.format(args.band, args.binsize)),
        np.transpose((psf.r[mask], psf.rprof[mask])))
p0 = (1e-2, 4, 1e-2, 6, 4)
coeff, var = curve_fit(_double_gauss, psf.r, psf.rprof, p0=p0)
#             sigma=err[mask])
pickle.dump( coeff, open( "psf_{0}.p".format(args.band), "wb" ) )
xfit = np.linspace(0, 40, 100)
yfit = _double_gauss(xfit, *coeff)
plt.plot(xfit, yfit, 'r')
np.savetxt(join(datadir, 'ascii', 'PSF_{0}_gauss.dat'.format(args.band)),
            np.transpose((xfit, yfit)))
plt.show()
