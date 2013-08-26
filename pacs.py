#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import glob
import pywcs
import argparse
from christensen import datadir, figsdir, horizons_file
import matplotlib.cm as cm
from scipy import ndimage
from datetime import datetime
from hsso import gildas

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--obsid', default=1342186621, type=int,
                choices=(1342186621, 1342203478))
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
parser.add_argument('-b', '--band', default="blue", choices=("blue", "red"),
                help="PACS band")
parser.add_argument('-s', '--save', default="store_true", help="save image")
parser.add_argument('--binsize', default=0, type=int, help="bin size")
args = parser.parse_args()

class Pacsmap(object):
    def __init__(self, obsid, zoom=0):
        """return patch centered on the nucleus"""
        self.fitsfile = glob.glob(join(datadir, str(obsid), 'level2',
            'HPPPMAP{}'.format(args.band[0].upper()), '*fits.gz'))[0]
        self.hdus = pyfits.open(self.fitsfile)
        self.cdelt2 = self.hdus[1].header['CDELT2']*3600

        # swap to little-endian byte order to avoid matplotlib bug
        pmap = self.hdus[1].data
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
        alpha = np.pi/2 + phase_ang*np.pi/180
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
        self.fov = int(round(60/self.pix))
        # patch with 2fovx2fov
        self.patch = pmap[com[1]-self.fov:com[1]+self.fov+1,
                          com[0]-self.fov:com[0]+self.fov+1]
        if zoom: self.patch = ndimage.zoom(self.patch, zoom, order=2)
        if args.debug:
            plt.imshow(pmap, origin="lower")
            plt.scatter(*comet)
            plt.show()
            plt.close()

    def shift(self):
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
        self.fov = int(round(30/self.pix))
        self.patch = self.patch[com[0]-self.fov:com[0]+self.fov+1,
                                com[1]-self.fov:com[1]+self.fov+1]

    def add(self, pmap):
        self.patch = np.average((self.patch, pmap.patch), axis=0)

    def radprof(self, binsize=0):
        """calculate radial profile"""
        y, x = np.indices(self.patch.shape)
        j, i = np.unravel_index(np.argmax(self.patch), self.patch.shape)
        r = np.sqrt((x-i)**2 + (y-j)**2)
        ind = np.argsort(r.flat)
        r *= self.cdelt2
        sr = r.flat[ind]
        sim = self.patch.flat[ind]
        if binsize:
            sr /= binsize
            ri = sr.astype(np.int16)
            deltar = ri[1:] - ri[:-1]
            rind = np.where(deltar)[0]
            rind = np.append([0], rind+1)
            error = [np.std(sim[lo:hi]) for lo,hi in zip(rind[:-1], rind[1:])]
            rprof = [np.mean(sim[lo:hi]) for lo,hi in zip(rind[:-1], rind[1:])]
            rad = binsize*(np.unique(ri)[:-1] + 0.5)
            return rad, rprof, error
        else:
            return sr, sim, np.zeros(len(sr))

# average orthogonal scans
pmap = Pacsmap(args.obsid)
pmap.add(Pacsmap(args.obsid+1))
pmap.shift()
patch = pmap.patch

plt.imshow(patch, origin="lower", interpolation=None)
plt.colorbar()
fov = patch.shape[0]/2
plt.title('{0} {1}'.format(args.obsid, args.band))
if args.band == "blue":
    levels = np.arange(-1.4, 0.1, 0.1) + .95*np.log10(np.abs(patch)).max()
else:
    levels = np.arange(-1., 0.1, 0.1) + .99*np.log10(np.abs(patch)).max()
zpatch = ndimage.zoom(patch, 4)
X = np.linspace(0, patch.shape[0]-1, len(zpatch))
plt.contour(X, X, np.log10(np.abs(zpatch)), levels=levels)
# plt.contour(np.log10(np.abs(patch)), levels=levels)
ax = plt.gca()
# ax.set_axis_off()
ax.autoscale(False)
plt.scatter(fov+pmap.comet[1], fov+pmap.comet[0], marker='x', color='k')
extent = ax.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
plt.savefig(join(figsdir, '{0}_{1}.png'.format(args.obsid, args.band)),
            bbox_inches=extent)
plt.show()
plt.close()

if args.binsize:
    r, prof, err = pmap.radprof(binsize=args.binsize)
    np.savetxt(join(datadir, 'ascii', '{0}_{1}_{2}_prof.dat'.format(args.obsid,
                args.band, args.binsize)),
                np.transpose((r, prof, err)))
    plt.errorbar(r, prof, yerr=err, fmt='o')
    for i in range(0, 30, args.binsize):
        plt.axvline(x=i, linestyle='--')
r, prof, err = pmap.radprof()
plt.scatter(r, prof, marker='x')
ax = plt.gca()
# ax.set_yscale('log')
# ax.set_xscale('log')
plt.xlim(0, 30)
# plt.ylim(1e-2, 0.4)
plt.show()
