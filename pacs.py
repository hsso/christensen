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
        pmap = self.hdus[1].data.byteswap().newbyteorder()
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
        pmap = ndimage.interpolation.shift(pmap, sh[::-1])
        pix = np.abs(self.cdelt2)
        fov = int(round(60/pix))
        self.patch = pmap[com[1]-fov:com[1]+fov+1, com[0]-fov:com[0]+fov+1]
        if zoom: self.patch = ndimage.zoom(self.patch, zoom, order=2)
        if args.debug:
            plt.imshow(pmap, origin="lower")
            plt.scatter(*comet)
            # plot max position
            mapmax = np.unravel_index(np.argmax(pmap), pmap.shape)[::-1]
            plt.scatter(*mapmax)
            # plot center-of-mass
            mapcom = ndimage.measurements.center_of_mass(pmap)[::-1]
            plt.scatter(*mapcom, color='r')
            plt.show()
            plt.close()

    def shift(self):
        plt.imshow(self.patch, origin="lower")
        mapmax = ndimage.measurements.maximum_position(self.patch)[::-1]
        plt.scatter(*mapmax)
        # plot center-of-mass
        mapcom = ndimage.measurements.maximum_position(self.patch)[::-1]
        plt.scatter(*mapcom, color='k')
        plt.show()
        plt.close()

    def add(self, pmap):
        self.patch = np.average((self.patch, pmap.patch), axis=0)

    def radprof(self, binsize=0):
        y, x = np.indices(self.patch.shape)
        j, i = np.unravel_index(np.argmax(self.patch), self.patch.shape)
        r = np.sqrt((x-i)**2 + (y-j)**2)
        ind = np.argsort(r.flat)
        r *= self.cdelt2
        sr = r.flat[ind]
        sim = self.patch.flat[ind]
        if binsize:
            ri = sr.astype(np.int16)
            deltar = ri[1:] - ri[:-1]
            rind = np.where(deltar)[0]
            nr = rind[1:] - rind[:-1]
            csim = np.cumsum(sim)
            tbin = csim[rind[1:]] - csim[rind[:-1]]
            rprof = tbin/nr
            return range(len(rprof)), rprof
        else:
            return sr, sim

# average orthogonal scans
pmap = Pacsmap(args.obsid)
pmap.add(Pacsmap(args.obsid+1))
patch = pmap.patch

if args.debug:
    plt.plot(patch.flat)
    plt.show()
plt.imshow(patch, origin="lower", cmap=cm.gist_heat_r)
plt.colorbar()
fov = patch.shape[0]/2
# plt.scatter(fov, fov)
plt.title('{0} {1}'.format(args.obsid, args.band))
if args.band == "blue":
    levels = np.arange(-1.4, 0.1, 0.1) + .99*np.log10(np.abs(patch)).max()
else:
    levels = np.arange(-1., 0.1, 0.1) + .99*np.log10(np.abs(patch)).max()
plt.contour(np.log10(np.abs(patch)), levels=levels, colors='g')
ax = plt.gca()
ax.set_axis_off()
extent = ax.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
plt.savefig(join(figsdir, '{0}_{1}.png'.format(args.obsid, args.band)),
            bbox_inches=extent)
plt.show()
plt.close()

r, prof = pmap.radprof(binsize=args.binsize)
np.savetxt(join(datadir, 'ascii', '{0}_{1}_{2}_prof.dat'.format(args.obsid,
            args.band, args.binsize)),
            np.transpose((r, prof)))
plt.scatter(r, prof)
ax = plt.gca()
# ax.set_yscale('log')
# ax.set_xscale('log')
plt.xlim(0, 30)
# plt.ylim(1e-2, 0.4)
plt.show()
