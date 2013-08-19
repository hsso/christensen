#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import glob
import pywcs
import argparse
from christensen import datadir, radec
import matplotlib.cm as cm

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--obsid', default=1342186621, type=int,
                choices=(1342186621, 1342186622, 1342203478, 1342203479))
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
args = parser.parse_args()

# pmap[pmap < cutlevels[i][0]] = cutlevels[i][0]

for c in ('Blue', 'Red'):
    fitsfile = glob.glob(join(datadir, str(args.obsid), 'level2',
        'HPPPMAP{}'.format(c[0]), '*fits.gz'))[0]
    hdulist = pyfits.open(fitsfile)

    # swap to little-endian byte order to avoid matplotlib bug
    pmap = hdulist[1].data.byteswap().newbyteorder()
    cdelt2 = hdulist[1].header['CDELT2']
    date_obs = hdulist[0].header['DATE-OBS']
    wcs = pywcs.WCS(hdulist[1].header)
    comet = wcs.wcs_sky2pix([radec], 1)[0]
    pix = np.abs(cdelt2)*3600
    fov = int(round(60/pix))
    print fov, date_obs
    com = [int(i) for i in comet]
#     patch = pmap[com[1]-fov:com[1]+fov, com[0]-fov:com[0]+fov]

    med = np.std(pmap[pmap > 0.])
    signal_stdev= np.std(pmap[pmap > med])
    print med, signal_stdev

#     plt.hist(pmap)
#     pmap[pmap > 0.1] = 0.1
#     pmap[pmap < 0.003] = -0.003
    if args.debug:
        plt.plot(pmap.flat)
        plt.show()
    plt.imshow(pmap, origin="lower", cmap=cm.gist_heat_r)
#     plt.scatter(*comet)
#     plt.scatter(fov, fov)
    plt.colorbar()
    plt.title('{0} {1}'.format(args.obsid, c))
#     plt.plot(ra, dec, 'ro-')
#     plt.contour(np.log10(patch), colors='k')
    plt.show()
