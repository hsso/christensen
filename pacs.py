#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import glob
import astropy.wcs as pywcs
import argparse
from christensen import datadir, figsdir, horizons_file
import matplotlib.cm as cm
from scipy import ndimage, stats, signal
from scipy.optimize import curve_fit
from datetime import datetime
from hsso import gildas
from herschel import Pacsmap
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--obsid', default=1342186621, type=int,
                choices=(1342186621, 1342203478))
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
parser.add_argument('-b', '--band', default="blue", choices=("blue", "red"),
                help="PACS band")
parser.add_argument('--profile', action="store_true", help="radial profile")
parser.add_argument('--save', action="store_true", help="save profile")
parser.add_argument('--fit', action="store_true", help="fit power law")
parser.add_argument('--binsize', default=0, type=float, help="bin size")
parser.add_argument('--rmax', default=0, type=int, help="radial profile extent")
args = parser.parse_args()

def mirror(arr, sign=1):
    return np.append(sign*arr[::-1], arr)

def _double_gauss(x, *p):
    """
    Sum of Gaussian functions
    """
    gauss1 = p[0]*stats.norm.pdf(x, 0, p[1])
    gauss2 = p[2]*stats.norm.pdf(x, p[3], p[4])
    return gauss1+gauss2

def _conv_prof(x, *p):
    """
    Convolved power law

    x : map profile
    """
    prof_con = signal.convolve(mirror(p[0]*pmap.r**-p[1] + p[2]), mirror(x),
                                mode="same")
    prof_con = ndimage.interpolation.shift(prof_con, -args.binsize*0.5,
                mode="nearest")
    if len(p)>3: prof_con += p[3]
    n = len(prof_con)/2
    if args.debug:
        plt.plot(mirror(pmap.r, sign=-1), prof_con)
        plt.plot(pmap.r, prof_con[n:])
    return prof_con[n:]

def pacsfile(obsid):
    fitsfile = glob.glob(join(datadir, str(obsid), 'level2',
        'HPPPMAP{}'.format(args.band[0].upper()), '*fits.gz'))[0]
    return fitsfile

if args.profile:
    # Use Pacsmap class in herschel module (untested)
    pmap = Pacsmap(pacsfile(args.obsid), fn=horizons_file[args.obsid],
            size=np.max((60, args.rmax+20)))
    pmap.add(Pacsmap(pacsfile(args.obsid+1), fn=horizons_file[args.obsid],
            size=np.max((60, args.rmax+20))))
    pmap.center_of_mass(size=40)
    pmap.shift(pmap.com, size=np.max((40, args.rmax)))
    fov = pmap.patch.shape[0]/2
    center = (fov+pmap.sh[1], fov+pmap.sh[0])
    if args.binsize:
        pmap.radprof(binsize=args.binsize, center=center, rmax=args.rmax)
        plt.errorbar(pmap.r, pmap.rprof, yerr=pmap.rprof_e, fmt='x', color="g")
        for i in np.arange(0, args.rmax, args.binsize):
            plt.axvline(x=i, linestyle='--')
    if args.fit:
        # read PSF gaussian coefficients
        coeff = pickle.load(open( "psf_{0}.p".format(args.band), "rb" ))
        p0 = (pmap.rprof[0], 1.1, 1e-3, 1e-4)
        popt, pcov = curve_fit(_conv_prof,
                                _double_gauss(pmap.r, *coeff),
                                pmap.rprof, p0=p0,
#                             sigma=pmap.rprof_e
                                )
        conv_prof = _conv_prof(
                            _double_gauss(pmap.r, *coeff),
                            *popt)
        print(popt, np.sqrt(pcov[1, 1]))
        print(conv_prof)
        plt.plot(pmap.r, conv_prof, color="green")
        np.savetxt(join(datadir, 'ascii',
                '{0}_{1}_{2}_{3}_prof_fit.dat'.format(args.obsid,
                args.band, args.binsize, args.rmax)),
                np.transpose((np.append(0, pmap.r), np.append(conv_prof[0],
                    conv_prof))))
    plt.scatter(mirror(pmap.r, sign=-1), mirror(pmap.rprof), marker='x',
            color='green')
    if args.save:
        np.savetxt(join(datadir, 'ascii',
                '{0}_{1}_{2}_{3}_prof.dat'.format(args.obsid,
                args.band, args.binsize, args.rmax)),
                np.transpose((pmap.r, pmap.rprof, pmap.rprof_e)))
    # unbinned profile
    pmap.radprof(center=center, rmax=args.rmax)
    plt.scatter(pmap.r, pmap.rprof, marker='.', color=args.band)
    if args.save:
        np.savetxt(join(datadir, 'ascii',
                '{0}_{1}_{2}_prof.dat'.format(args.obsid,
                args.band, args.rmax)),
                np.transpose((pmap.r, pmap.rprof, pmap.rprof_e)))
    ax = plt.gca()
    ax.set_yscale('log')
    if args.rmax: plt.xlim(0, args.rmax)
    plt.ylim(ymin=1e-6)
    plt.show()
else:
    pmap = Pacsmap(pacsfile(args.obsid), fn=horizons_file[args.obsid])
    pmap.add(Pacsmap(pacsfile(args.obsid+1), fn=horizons_file[args.obsid]))
    pmap.center_of_mass(size=30)
    pmap.shift(pmap.com, size=30)
    patch = pmap.patch
    plt.imshow(patch, origin="lower", interpolation=None, cmap=cm.gist_heat_r)
    plt.colorbar()
    fov = patch.shape[0]/2
    plt.title('{0} {1}'.format(args.obsid, args.band))
    # plot contours
    if args.band == "blue":
        levels = np.arange(-1.3, 0.1, 0.1) + .95*np.log10(np.abs(patch)).max()
    else:
        levels = np.arange(-0.7, 0.1, 0.1) + .95*np.log10(np.abs(patch)).max()
    zpatch = ndimage.zoom(patch, 4)
    X = np.linspace(0, patch.shape[0]-1, len(zpatch))
    plt.contour(X, X, np.log10(np.abs(zpatch)), levels=levels, zorder=1)
    # plt.contour(np.log10(np.abs(patch)), levels=levels)
    ax = plt.gca()
    ax.set_axis_off()
    ax.autoscale(False)
    plt.scatter(fov+pmap.comet[1], fov+pmap.comet[0], marker='o', color='y')
    extent = ax.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
    plt.savefig(join(figsdir, '{0}_{1}.png'.format(args.obsid, args.band)),
                bbox_inches=extent)
    plt.show()
    pmap.tofits(join(datadir, '{0}_{1}.fits'.format(args.obsid, args.band)))
