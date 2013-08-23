#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import glob
import pywcs
import argparse
from christensen import datadir, figsdir
import matplotlib.cm as cm
from scipy.ndimage import interpolation, zoom
from datetime import datetime
from hsso import gildas

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--obsid', default=1342186621, type=int,
                choices=(1342186621, 1342186622, 1342203478, 1342203479))
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
parser.add_argument('-b', '--band', default="blue", choices=("blue", "red"),
                help="PACS band")
parser.add_argument('-s', '--save', default="store_true", help="save image")
args = parser.parse_args()

# pmap[pmap < cutlevels[i][0]] = cutlevels[i][0]

def pacsmap(obsid):
    fitsfile = glob.glob(join(datadir, str(obsid), 'level2',
        'HPPPMAP{}'.format(args.band[0].upper()), '*fits.gz'))[0]
    hdus = pyfits.open(fitsfile)

    cdelt2 = hdus[1].header['CDELT2']
    date_obs = hdus[0].header['DATE-OBS']
    date_end = hdus[0].header['DATE-END']
    # calculate comet position at midtime
    start = datetime.strptime(date_obs[:20], "%Y-%m-%dT%H:%M:%S.")
    end = datetime.strptime(date_end[:20], "%Y-%m-%dT%H:%M:%S.")
    mid_time = start + (end-start)/2
    # interpolate ra and dec of the comet
    ra = gildas.deltadot(mid_time, filename="horizons_2009.txt", column=2)
    dec = gildas.deltadot(mid_time, filename="horizons_2009.txt", column=3)
    wcs = pywcs.WCS(hdus[1].header)
    # origin coordinate is 0 (Numpy and C standards)
    comet = wcs.wcs_sky2pix([(ra, dec)], 0)[0]
    # swap to little-endian byte order to avoid matplotlib bug
    pmap = hdus[1].data.byteswap().newbyteorder()
    com = [int(round(i)) for i in comet]
    sh  = [i - round(i) for i in comet]
    pmap = interpolation.shift(pmap, sh)
    pix = np.abs(cdelt2)*3600
    fov = int(round(30/pix))
    print comet, com, sh, fov, date_obs
    patch = pmap[com[1]-fov:com[1]+fov+1, com[0]-fov:com[0]+fov+1]
    return zoom(patch, 3, order=2)

# average orthogonal scans
pmap = np.average((pacsmap(args.obsid), pacsmap(args.obsid+1)), axis=0)

if args.debug:
    plt.plot(pmap.flat)
    plt.show()
plt.imshow(pmap, origin="lower", cmap=cm.gist_heat_r)
plt.colorbar()
# plt.scatter(*comet)
fov = pmap.shape[0]/2
# plt.scatter(fov, fov)
plt.title('{0} {1}'.format(args.obsid, args.band))
if args.band == "blue":
    levels = np.arange(-1.1, 0.1, 0.1) + .99*np.log10(np.abs(pmap)).max()
else:
    levels = np.arange(-.7, 0.1, 0.1) + .99*np.log10(np.abs(pmap)).max()
plt.contour(np.log10(np.abs(pmap)), levels=levels, colors='g')
ax = plt.gca()
ax.set_axis_off()
extent = ax.get_window_extent().transformed(plt.gcf().dpi_scale_trans.inverted())
plt.savefig(join(figsdir, '{0}_{1}.png'.format(args.obsid, args.band)),
            bbox_inches=extent)
plt.show()
