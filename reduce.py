#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser, join
import glob
import argparse
from scipy import interpolate
from hsso import gildas

def ff(pol):
    hdulist = pyfits.open( glob.glob(
            join(datadir, str(obsid), 'level2',
            '{0}-{1}-{2}'.format(args.backend, pol, args.sideband),
            'box_001', '*.fits*'))[0])

    for i in hdulist[1].header.ascardlist().keys():
        if hdulist[1].header[i] == 'loThrow':
            throw = hdulist[1].header[i[4:]]
            break
    freq = hdulist[1].data.field('{0}frequency_{1}'.format(args.sideband.lower(), args.subband))[0]
    flux = hdulist[1].data.field('flux_{0}'.format(args.subband))[0]
    return freq, flux, throw

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--backend', default='WBS', choices=('HRS', 'WBS'))
parser.add_argument('--sideband', default='USB', choices=('LSB', 'USB'))
parser.add_argument('--subband', default=1, type=int, choices=range(1,5))
args = parser.parse_args()

datadir = expanduser('~/HssO/Christensen/data')
obsid = 1342204014
freq0 = 556.9359877

freqh, fluxh, throwh = ff('H')
freqv, fluxv, throwv = ff('V')

freq_list = (freqh, freqh+throwh, freqv, freqv+throwv)
flux_list = (fluxh, -fluxh, fluxv, -fluxv)

freqav, fluxav = gildas.averagen(freq_list, flux_list)

plt.plot(freqav, fluxav,  drawstyle='steps-mid')
plt.axvline(x=freq0, linestyle='--')
plt.show()
