#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser, join
import glob
import argparse
from scipy import fftpack
from hsso import gildas
from herschel import fft

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

freq_list, flux_list = [], []
for pol in ('H', 'V'):
    hdulist = pyfits.open( glob.glob(
        join(datadir, str(obsid), 'level2',
        '{0}-{1}-{2}'.format(args.backend, pol, args.sideband),
        'box_001', '*.fits*'))[0])
    freq, flux, throw = fft(hdulist, args.sideband, args.subband)
    freq_list.extend([freq, freq+throw])
    flux_list.extend([flux, -flux])

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
