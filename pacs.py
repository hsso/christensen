#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import glob
import pywcsgrid2
import argparse
from christensen import datadir
import matplotlib.cm as cm

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--obsid', default=1342186621, type=int)
args = parser.parse_args()

for c in ('Blue', 'Red'):
    fitsfile = glob.glob(join(datadir, str(args.obsid), 'level2',
        'HPPPMAP{}'.format(c[0]), '*fits.gz'))[0]
    hdulist = pyfits.open(fitsfile)

    pmap = hdulist[1].data.byteswap().newbyteorder()

    pywcsgrid2.subplot(111, header=hdulist[1].header)
    plt.imshow(pmap, origin="lower", cmap=cm.gist_heat_r)
    plt.colorbar()
    plt.show()
