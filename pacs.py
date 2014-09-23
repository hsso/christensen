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

class Pacsmap2(object):
    """
    Read PACS photometry map
    """

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
            start = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S.%f")
            end = datetime.strptime(date_end, "%Y-%m-%dT%H:%M:%S.%f")
            mid_time = start + (end-start)/2
            # interpolate ra and dec of the comet
            ra = gildas.deltadot(mid_time, filename=horizons_file[args.obsid],
                    column=2)
            dec = gildas.deltadot(mid_time, filename=horizons_file[args.obsid],
                    column=3)
            # calculate direction toward the Sun
            phase_ang = gildas.deltadot(mid_time,
                    filename=horizons_file[args.obsid], column=8)
            alpha = 3*np.pi/2 - phase_ang*np.pi/180
            cos, sin = np.cos(alpha), np.sin(alpha)
            # origin coordinate is 0 (Numpy and C standards)
            wcs = pywcs.WCS(self.hdus[1].header)
            comet = wcs.wcs_world2pix([(ra, dec)], 0)[0]
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
#             plt.scatter(*comet)
            plt.show()
            plt.close()

    def com(self, size=30):
        """select the top 0.99% pixels to calculate the center of mass"""
        fov = int(round(size/self.pix))
        center = self.patch.shape[0]/2
        zoom = self.patch[center-fov:center+fov+1,
                            center-fov:center+fov+1]
        hist, bins = np.histogram(zoom.ravel(), normed=True, bins=100)
        threshold = bins[np.cumsum(hist) * (bins[1] - bins[0]) > 0.992][0]
        mpatch = np.ma.masked_less(zoom, threshold)
        mapcom = ndimage.measurements.center_of_mass(mpatch)
        if args.debug:
            plt.imshow(self.patch, origin="lower")
            mapmax = ndimage.measurements.maximum_position(zoom)[::-1]
            plt.scatter(*mapmax)
            # plot center-of-mass
            plt.scatter(*mapcom[::-1], color='r')
            plt.show()
            plt.close()
        # center of mass
        self.com = [int(round(i)) for i in mapcom]
        # fraction of pixel to com
        self.sh  = np.array(mapcom) - self.com
        self.com = [center - fov + i for i in mapcom]

    def shift(self, center, size=30):
        """shift array to be centerd at cneter"""
        self.comet = [self.fov-center[0], self.fov-center[1]]
        self.fov = int(round(size/self.pix))
        self.patch = self.patch[center[0]-self.fov:center[0]+self.fov+1,
                                center[1]-self.fov:center[1]+self.fov+1]

    def gauss_fit(self):
        """fit 2D Gaussian"""
        pass

    def add(self, pmap):
        """average orthogonal scans"""
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
        # normalize to Jy arcsec-2
        sim /= self.cdelt2**2
        if binsize:
            sr /= binsize
            ri = sr.astype(np.int16)
            deltar = ri[1:] - ri[:-1]
            rind = np.where(deltar)[0]
            rind = np.concatenate(([0], rind+1, [len(ri)]))
            n = rind[1:] - rind[:-1]
            self.rprof = np.array([np.mean(sim[lo:hi]) for lo,hi in
                                    zip(rind[:-1], rind[1:])])
            self.rprof_e = np.array([np.std(sim[lo:hi]) for lo,hi in
                                    zip(rind[:-1], rind[1:])])
            self.rprof_e = np.where(self.rprof > self.rprof_e, self.rprof_e,
                                    0.99*self.rprof)
            self.r = binsize*(np.unique(ri)[:-1] + 0.5)
        else:
            self.r = sr
            self.rprof = sim
            self.rprof_e = np.zeros(len(sr))
        if rmax:
            mask = [self.r < rmax]
            self.r = self.r[mask]
            self.rprof = self.rprof[mask]
            self.rprof_e = self.rprof_e[mask]

if args.profile:
    pmap = Pacsmap(args.obsid, size=np.max((60, args.rmax+20)))
    pmap.add(Pacsmap(args.obsid+1, size=np.max((60, args.rmax+20))))
    pmap.com(size=40)
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
    pmap.com(size=30)
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
