#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import glob
import pywcs
import argparse
from christensen import datadir, figsdir, horizons_file, psfdir, psf_vesta
import matplotlib.cm as cm
from scipy import ndimage, stats, signal
from scipy.optimize import curve_fit
from datetime import datetime
from hsso import gildas

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--obsid', default=1342186621, type=int,
                choices=(1342186621, 1342203478))
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
parser.add_argument('-b', '--band', default="blue", choices=("blue", "red"),
                help="PACS band")
parser.add_argument('--profile', action="store_true", help="radial profile")
parser.add_argument('--psf', action="store_true", help="PSF profile")
parser.add_argument('--save', action="store_true", help="save profile")
parser.add_argument('--binsize', default=0, type=float, help="bin size")
parser.add_argument('--rmax', default=0, type=int, help="radial profile extent")
args = parser.parse_args()

def _double_gauss(x, *p):
    """Gaussian fitting"""
    gauss1 = p[0]*stats.norm.pdf(x, 0, p[1])
    gauss2 = p[2]*stats.norm.pdf(x, p[3], p[4])
    return gauss1+gauss2

class Pacsmap(object):
    def __init__(self, obsid, size=60, zoom=0, comet=True):
        """return patch centered on the nucleus"""
        if isinstance(obsid, int):
            self.fitsfile = glob.glob(join(datadir, str(obsid), 'level2',
                'HPPPMAP{}'.format(args.band[0].upper()), '*fits.gz'))[0]
        else:
            self.fitsfile = obsid
        self.hdus = pyfits.open(self.fitsfile)
        self.cdelt2 = self.hdus[1].header['CDELT2']*3600

        # swap to little-endian byte order to avoid matplotlib bug
        pmap = self.hdus[1].data
        if comet:
            # calculate comet position at midtime
            date_obs = self.hdus[0].header['DATE-OBS']
            date_end = self.hdus[0].header['DATE-END']
            start = datetime.strptime(date_obs[:20], "%Y-%m-%dT%H:%M:%S.")
            end = datetime.strptime(date_end[:20], "%Y-%m-%dT%H:%M:%S.")
            mid_time = start + (end-start)/2
            # interpolate ra and dec of the comet
            ra = gildas.deltadot(mid_time, filename=horizons_file[args.obsid], column=2)
            dec = gildas.deltadot(mid_time, filename=horizons_file[args.obsid], column=3)
            # calculate direction toward the Sun
            phase_ang = gildas.deltadot(mid_time, filename=horizons_file[args.obsid], column=8)
            alpha = 3*np.pi/2 - phase_ang*np.pi/180
            cos, sin = np.cos(alpha), np.sin(alpha)
            print("\draw[->,yellow] (axis cs:{0:.3f},{1:.3f}) --\n"
                    "(axis cs:{2:.3f},{3:.3f});".format(10*cos, 10*sin, 20*cos, 20*sin))
            print(r"\node[yellow] at (axis cs:{0:.3f},{1:.3f}) {{\sun}};".format(25*cos,
                            25*sin))
            # origin coordinate is 0 (Numpy and C standards)
            wcs = pywcs.WCS(self.hdus[1].header)
            comet = wcs.wcs_sky2pix([(ra, dec)], 0)[0]
            com = [int(round(i)) for i in comet]
            sh  = comet-com
            # shift array to center on comet nucleus
            pmap = ndimage.interpolation.shift(pmap, sh[::-1])
            self.pix = np.abs(self.cdelt2)
            self.fov = int(round(size/self.pix))
            # patch with 2fovx2fov
            self.patch = pmap[com[1]-self.fov:com[1]+self.fov+1,
                              com[0]-self.fov:com[0]+self.fov+1]
        else:
            self.patch = pmap
        if zoom: self.patch = ndimage.zoom(self.patch, zoom, order=2)
        if args.debug:
            plt.imshow(pmap, origin="lower")
            plt.scatter(*comet)
            plt.show()
            plt.close()

    def shift(self, size=30):
        # select the top 0.99% pixels to calculate the center of mass
        hist, bins = np.histogram(self.patch.ravel(), normed=True, bins=100)
        threshold = bins[np.cumsum(hist) * (bins[1] - bins[0]) > 0.992][0]
        mpatch = np.ma.masked_less(self.patch, threshold)
        mapcom = ndimage.measurements.center_of_mass(mpatch)
        if args.debug:
            plt.imshow(self.patch, origin="lower")
            mapmax = ndimage.measurements.maximum_position(self.patch)[::-1]
            plt.scatter(*mapmax)
            # plot center-of-mass
            plt.scatter(*mapcom[::-1], color='r')
            plt.show()
            plt.close()
        com = [int(round(i)) for i in mapcom]
        self.comet = [self.fov-com[0], self.fov-com[1]]
        self.sh  = np.array(mapcom)-com
#         self.patch = ndimage.interpolation.shift(self.patch, sh)
        self.fov = int(round(size/self.pix))
        self.patch = self.patch[com[0]-self.fov:com[0]+self.fov+1,
                                com[1]-self.fov:com[1]+self.fov+1]

    def add(self, pmap):
        self.patch = np.average((self.patch, pmap.patch), axis=0)

    def radprof(self, center=None, binsize=0, rmax=0):
        """calculate radial profile"""
        y, x = np.indices(self.patch.shape)
        if center:
            i, j = center
        else:
            j, i = np.unravel_index(np.argmax(self.patch), self.patch.shape)
        r = np.sqrt((x-i)**2 + (y-j)**2)
        ind = np.argsort(r.flat)
        r *= self.cdelt2
        sr = r.flat[ind]
        sim = self.patch.flat[ind]
        sr = sr[sim >0]
        sim = sim[sim >0]
        if binsize:
            sr /= binsize
            ri = sr.astype(np.int16)
            deltar = ri[1:] - ri[:-1]
            rind = np.where(deltar)[0]
            rind = np.append([0], rind+1)
            self.rprof = np.array([np.mean(sim[lo:hi]) for lo,hi in
                                    zip(rind[:-1], rind[1:])])
            self.rprof_e = np.array([np.std(sim[lo:hi]) for lo,hi in
                                    zip(rind[:-1], rind[1:])])
#             self.rprof_e = np.where(self.rprof > self.rprof_e, self.rprof_e, 0.99*self.rprof)
            self.rad = binsize*(np.unique(ri)[:-1] + 0.5)
        else:
            self.rad = sr
            self.rprof = sim
            self.rprof_e = np.zeros(len(sr))
        if rmax:
            mask = [self.rad < rmax]
            self.rad = self.rad[mask]
            self.rprof = self.rprof[mask]
            self.rprof_e = self.rprof_e[mask]

def mirror(arr, sign=1):
    return np.append(sign*arr[::-1], arr)

psf = Pacsmap(join(psfdir, psf_vesta[args.band]), comet=False,
                size=np.max((60, args.rmax)))
fov = psf.patch.shape[0]/2
psf.radprof(center=(fov,fov), rmax=args.rmax)
plt.scatter(psf.rad, psf.rprof, marker='x', color=args.band)
ax = plt.gca()
ax.set_yscale('log')
if args.debug:
    plt.imshow(psf.patch, origin="lower")
    plt.scatter(fov,fov)
    plt.show()
if args.binsize:
    psf.radprof(binsize=args.binsize, center=(fov,fov), rmax=args.rmax)
    plt.errorbar(psf.rad, psf.rprof, yerr=psf.rprof_e, fmt='x', color="g")
    for i in np.arange(0, args.rmax, args.binsize):
        plt.axvline(x=i, linestyle='--')
    np.savetxt(join(datadir, 'ascii', 'PSF_{0}_{1}_psfprof.dat'.format(args.band,
            args.binsize)),
            np.transpose((psf.rad, psf.rprof, psf.rprof_e)))
    plt.scatter(mirror(psf.rad, sign=-1), mirror(psf.rprof), marker='x', color="g")
# plt.xlim(0,40)
# plt.show()

if args.profile:
    # average orthogonal scans
    pmap = Pacsmap(args.obsid, size=np.max((60, args.rmax+20)))
    pmap.add(Pacsmap(args.obsid+1, size=np.max((60, args.rmax+20))))
    pmap.shift(size=np.max((40, args.rmax)))
    patch = pmap.patch
    fov = patch.shape[0]/2
    center = (fov+pmap.sh[1], fov+pmap.sh[0])
    if args.binsize:
        pmap.radprof(binsize=args.binsize, center=center, rmax=args.rmax)
        if args.save:
            np.savetxt(join(datadir, 'ascii', '{0}_{1}_{2}_prof.dat'.format(args.obsid,
                    args.band, args.binsize)),
                    np.transpose((pmap.rad, pmap.rprof, pmap.rprof_e)))
        plt.errorbar(pmap.rad, pmap.rprof, yerr=pmap.rprof_e, fmt='x', color="g")
        for i in np.arange(0, args.rmax, args.binsize):
            plt.axvline(x=i, linestyle='--')
#     prof_decon, error = signal.deconvolve(mirror(pmap.rprof), mirror(psf.rprof))
    prof_con = signal.convolve(mirror(pmap.rad)**-1.05, mirror(psf.rprof),
                                mode="same")
    plt.plot(mirror(pmap.rad, sign=-1)-args.binsize*.5,
                    pmap.rprof[0]*prof_con, color="red")
    plt.scatter(mirror(pmap.rad, sign=-1), mirror(pmap.rprof), marker='x', color='green')
    pmap.radprof(center=center, rmax=args.rmax)
    plt.scatter(pmap.rad, pmap.rprof, marker='x', color=args.band)
    if args.save:
        np.savetxt(join(datadir, 'ascii', '{0}_{1}_prof.dat'.format(args.obsid,
                args.band)),
                np.transpose((pmap.rad, pmap.rprof, pmap.rprof_e)))
    ax = plt.gca()
#     ax.set_yscale('log')
    if args.rmax: plt.xlim(-args.rmax, args.rmax)
#     plt.ylim(1e-5)
    plt.show()
else:
    # average orthogonal scans
    pmap = Pacsmap(args.obsid)
    pmap.add(Pacsmap(args.obsid+1))
    pmap.shift(size=30)
    patch = pmap.patch
    plt.imshow(patch, origin="lower", interpolation=None, cmap=cm.gist_heat_r)
    plt.colorbar()
    fov = patch.shape[0]/2
    plt.title('{0} {1}'.format(args.obsid, args.band))
    if args.band == "blue":
        levels = np.arange(-1.3, 0.1, 0.1) + .95*np.log10(np.abs(patch)).max()
    else:
        levels = np.arange(-0.8, 0.1, 0.1) + .95*np.log10(np.abs(patch)).max()
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
