#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser, join
import glob
import argparse
from scipy import interpolate, fftpack
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
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
parser.add_argument('--fftlim', default=2e2, type=float,
                    help='FFT high frequency limit')
args = parser.parse_args()

datadir = expanduser('~/HssO/Christensen/data')
obsid = 1342204014
freq0 = 556.9359877

freqh, fluxh, throwh = ff('H')
freqv, fluxv, throwv = ff('V')

freq_list = (freqh, freqh+throwh, freqv, freqv+throwv)
flux_list = (fluxh, -fluxh, fluxv, -fluxv)

freqav, fluxav = gildas.averagen(freq_list, flux_list, goodval=True)
vel = gildas.vel(freqav, freq0)

sample_freq = fftpack.fftfreq(fluxav.size, d=np.abs(freqav[0]-freqav[1]))
sig_fft = fftpack.fft(fluxav)

if args.debug:
    pidxs = np.where(sample_freq > 0)
    f = sample_freq[pidxs]
    pgram = np.abs(sig_fft)[pidxs]
    plt.loglog(f, pgram)
    plt.axvline(x=args.fftlim, linestyle='--')
    plt.show()

sig_fft[np.abs(sample_freq) > args.fftlim] = 0
baseline = np.real(fftpack.ifft(sig_fft))

if args.debug:
    plt.plot(freqav, fluxav,  drawstyle='steps-mid')
    plt.plot(freqav, baseline)
    plt.axvline(x=freq0, linestyle='--')
    plt.show()

fluxav -= baseline
print('rms = {0:.2f} mK'.format(np.std(fluxav[4:-4])*1e3))
mask = [np.abs(vel) <= 20]
plt.plot(vel, fluxav,  drawstyle='steps-mid')
plt.show()

np.savetxt(expanduser("~/HssO/Christensen/data/ascii/{}_{:.0f}_{}.dat".format(
                    obsid, freq0, args.backend)),
                    np.transpose((vel[mask], fluxav[mask])))
