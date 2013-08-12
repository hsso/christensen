#!/usr/bin/python

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser, join
import glob
import argparse
from scipy import linalg, fftpack
from hsso import gildas
from hsso.class_utils import pgram_peaks, linfunc, fitfunc
from herschel import fft
from christensen import datadir, freq0

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--backend', default='WBS', choices=('HRS', 'WBS'))
parser.add_argument('--sideband', default='', choices=('', 'LSB', 'USB'))
parser.add_argument('--subband', default=0, type=int, choices=range(5))
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
parser.add_argument('--fftlim', default=2e2, type=float,
                    help='FFT high frequency limit')
parser.add_argument('-n', '--num', default=8, type=int,
                    help='number of sine waves')
parser.add_argument('-m', '--mol', default='H2O', choices=('H2O', 'NH3'))
args = parser.parse_args()

subband = {'HRS': 1, 'WBS': 4}
sideband = {'H2O': 'LSB', 'NH3': 'USB'}
if not args.subband: args.subband = subband[args.backend]
if not args.sideband: args.sideband = sideband[args.mol]

obsid = 1342204014

freq_list, flux_list = [], []
for pol in ('H', 'V'):
    """Return list of frequencies and fluxes in little endian byte order"""
    hdulist = pyfits.open( glob.glob(
        join(datadir, str(obsid), 'level2',
        '{0}-{1}-{2}'.format(args.backend, pol, args.sideband),
        'box_001', '*.fits*'))[0])
    freq, flux, throw = fft(hdulist, args.sideband, args.subband)
    flux = flux.byteswap().newbyteorder('L')/.75
    freq = freq.byteswap().newbyteorder('L')
    freq_list.extend([freq, freq+throw])
    flux_list.extend([flux, -flux])

freqav, fluxav = gildas.averagen(freq_list, flux_list, goodval=True)
vel = gildas.vel(freqav, freq0[args.mol])

# FFT
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
    plt.axvline(x=freq0[args.mol], linestyle='--')
    plt.show()

fluxav -= baseline

# Lomb-Scargle periodogram
scaled_flux = fluxav-fluxav.mean()
# frequency
f = np.linspace(1e2, 2e4, 1e4)
pgram, peak_freqs, peak_flux = pgram_peaks(freqav, scaled_flux, f, args.num)
if args.debug:
    plt.loglog(f, pgram)
    for maxfreq in peak_freqs:
        plt.axvline(x=maxfreq, linestyle='--')
    plt.show()

A = linfunc(np.ones(2*args.num), freqav, peak_freqs)
c, resid, rank, sigma = linalg.lstsq(A, scaled_flux)
baseline = np.sum(A*c, axis=1) + fluxav.mean()

if args.debug:
    plt.plot(freqav, fluxav,  drawstyle='steps-mid')
    plt.plot(freqav, baseline)
    plt.axvline(x=freq0[args.mol], linestyle='--')
#     plt.plot(freqav, 0.03*np.sin(np.max(peak_freqs)*freqav), 'red')
    plt.show()

fluxav -= baseline
delv = np.abs(np.average(vel[1:]-vel[:-1]))
n = np.ceil(2*0.4/delv)
rms = np.std(fluxav[4:-4])*1e3
upper = 3 * np.sqrt(n) * delv * rms
print('rms = {0:.2f} mK, delv = {1:.3f} km/s, upper= {2:.2f} K m/s'.format(rms,
        delv, upper))
mask = [np.abs(vel) <= 20]
plt.plot(vel, fluxav,  drawstyle='steps-mid')
plt.show()

np.savetxt(expanduser("~/HssO/Christensen/data/ascii/{}_{:.0f}_{}.dat".format(
                    obsid, freq0[args.mol], args.backend)),
                    np.transpose((vel[mask], fluxav[mask]*1e3)))
