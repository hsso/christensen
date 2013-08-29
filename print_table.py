#!/usr/bin/python

from christensen import obsids, datadir, paperdir
from os.path import join
import glob
import pyfits

f = open(join(paperdir, 'log.txt'), 'w')
for obsid in obsids:
    hdus = pyfits.open(glob.glob(join(datadir, str(obsid), '*fits.gz'))[0])
    print hdus[0].header["DATE-OBS"]
    f.write("& {0} \\\\\n".format(obsid))
