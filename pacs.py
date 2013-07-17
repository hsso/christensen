#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import glob
import pywcsgrid2
import pywcs
import argparse
from christensen import *
import matplotlib.cm as cm

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--obsid', default=1342186621, type=int)
args = parser.parse_args()

for c in ('Blue', 'Red'):
    fitsfile = glob.glob(join(datadir, str(args.obsid), 'level2',
        'HPPPMAP{}'.format(c[0]), '*fits.gz'))[0]
    hdulist = pyfits.open(fitsfile)

    pmap = hdulist[1].data.byteswap().newbyteorder()
    sh = pmap.shape
    patch = pmap[sh[0]/2-5:sh[0]/2+5, sh[1]/2-5:sh[1]/2+5]
    h = hdulist[1].header
    wcs = pywcs.WCS(hdulist[1].header)
    comet = wcs.wcs_sky2pix([radec1, radec2], 1)
    print hdulist[0].header['DATE-OBS']

    ax = pywcsgrid2.subplot(111, header= h)
    ax.set_ticklabel_type("delta", "delta", offset=radec[0],
                    latitude=radec[1])
#     ax.grid()
    plt.imshow(pmap, origin="lower", cmap=cm.gist_heat_r)
#     plt.plot(ra, dec, 'ro-')
    for i in comet: plt.scatter(*i)
#     plt.colorbar()
#     plt.contour(np.log10(patch), colors='k')
    plt.show()
