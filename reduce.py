#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser, join
import glob
import argparse
from scipy import linalg, fftpack
from hsso import gildas
from hsso.class_utils import pgram_peaks, linfunc, fitfunc
from herschel import HIFISpectrum
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
parser.add_argument('--deg', default=2, type=int,
                    help='baseline polynomial degree')
parser.add_argument('--twiny', action='store_true')
parser.add_argument("--lim", nargs=2, type=float, default=(-1, 1),
                    help='line limits in km/s')
args = parser.parse_args()

subband = {'HRS': 1, 'WBS': 4}
sideband = {'H2O': 'LSB', 'NH3': 'USB'}
if not args.subband: args.subband = subband[args.backend]
if not args.sideband: args.sideband = sideband[args.mol]

obsid = 1342204014

freq_list, flux_list = [], []
for pol in ('H', 'V'):
    """Return list of frequencies and fluxes in little endian byte order"""
    spec = HIFISpectrum(glob.glob(
        join(datadir, str(obsid), 'level2',
        '{0}-{1}-{2}'.format(args.backend, pol, args.sideband),
        'box_001', '*.fits*'))[0], args.sideband, args.subband)
    freq_list.extend([spec.freq, spec.freq+spec.throw])
    flux_list.extend([spec.flux, -spec.flux])
    spec.save(join(datadir, 'ascii'))
    print gildas.vel(freq0[args.mol] - spec.throw, freq0[args.mol])
    plt.plot(spec.freq, spec.flux, drawstyle='steps-mid')
    plt.axvline(x=freq0[args.mol], linestyle='--')
    plt.axvline(x=freq0[args.mol] - spec.throw, linestyle='--')
    plt.show()

freqav, fluxav = gildas.averagen(freq_list, flux_list, goodval=True)
vel = gildas.vel(freqav, freq0[args.mol])

baseflux = fluxav.copy()
maskline = np.where(np.abs(vel) < 1)
maskvel = np.where((np.abs(vel) < 3) & (np.abs(vel) > 1))
velmask = np.where((np.abs(vel) < 3))
func = np.poly1d(np.polyfit(freqav[maskvel], baseflux[maskvel], args.deg))
baseflux[maskline] = func(freqav[maskline])

# FFT
sample_freq = fftpack.fftfreq(fluxav.size, d=np.abs(freqav[0]-freqav[1]))
sig_fft = fftpack.fft(baseflux)

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
    plt.plot(freqav[velmask], func(freqav[velmask]))
    plt.axvline(x=freq0[args.mol], linestyle='--')
    if args.twiny:
        ax1 = plt.gca()
        # update xlim of ax2
        ax2 = ax1.twiny()
        x1, x2 = ax1.get_xlim()
        ax2.set_xlim(gildas.vel(x1, freq0[args.mol]),
                     gildas.vel(x2, freq0[args.mol]))
    plt.show()

baseflux -= baseline

# Lomb-Scargle periodogram
scaled_flux = baseflux-baseflux.mean()
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
lomb_baseline = np.sum(A*c, axis=1) + baseflux.mean()
baseline += lomb_baseline

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
intens = gildas.intens(fluxav, vel, lim=args.lim)
print('intens = {0[0]} {0[1]} K km/s'.format(intens))
print('rms = {0:.2f} mK, delv = {1:.3f} km/s, upper= {2:.2f} K m/s'.format(rms,
        delv, upper))
mask = [np.abs(vel) <= 20]
plt.plot(vel[mask], fluxav[mask],  drawstyle='steps-mid')
# plt.axvspan(*args.lim, facecolor='b', alpha=0.5)
plt.axhline(y=0)
plt.axvline(x=0)
plt.savefig('{0}_{1}_baseline.pdf'.format(obsid, args.backend))
plt.show()

np.savetxt(expanduser("~/HssO/Christensen/data/ascii/{}_{:.0f}_{}.dat".format(
                    obsid, freq0[args.mol], args.backend)),
                    np.transpose((vel[mask], fluxav[mask]*1e3)))
