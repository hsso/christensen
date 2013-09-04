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

def fitsfile(obsid, backend, pol, sideband):
    return glob.glob(
        join(datadir, str(obsid), 'level2',
        '{0}-{1}-{2}'.format(backend, pol, sideband),
        'box_001', '*.fits*'))[0]

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

spec = HIFISpectrum(fitsfile(obsid, args.backend, 'H', args.sideband),
                    args.subband)
spec.scale((-60, 10))
spec.save(join(datadir, 'ascii'))
specv = HIFISpectrum(fitsfile(obsid, args.backend, 'V', args.sideband),
                    args.subband)
specv.scale((-60, 10))
specv.save(join(datadir, 'ascii'))
# spec_ave = spec + specv
# spec_ave.save(join(datadir, 'ascii'), '_aver')
spec.fold()
spec.scale((-10, 10))
spec.save(join(datadir, 'ascii'), '_folded')
specv.fold()
specv.scale((-10, 10))
specv.save(join(datadir, 'ascii'), '_folded')
spec.add(specv)
spec.scale((-10, 10))
spec.save(join(datadir, 'ascii'), '_ave')
spec.baseline(args.fftlim, join(datadir, 'ascii'))

baseflux = spec.flux.copy()
maskline = np.where(np.abs(spec.vel) < 1)
maskvel = np.where((np.abs(spec.vel) < 3) & (np.abs(spec.vel) > 1))
velmask = np.where((np.abs(spec.vel) < 3))
func = np.poly1d(np.polyfit(spec.freq[maskvel], baseflux[maskvel], args.deg))
baseflux[maskline] = func(spec.freq[maskline])

# FFT
sample_freq = fftpack.fftfreq(spec.flux.size, d=np.abs(spec.freq[0]-spec.freq[1]))
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
    plt.plot(spec.freq, spec.flux,  drawstyle='steps-mid')
    plt.plot(spec.freq, baseline)
    plt.plot(spec.freq[velmask], func(spec.freq[velmask]))
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
pgram, peak_freqs, peak_flux = pgram_peaks(spec.freq, scaled_flux, f, args.num)
if args.debug:
    plt.loglog(f, pgram)
    for maxfreq in peak_freqs:
        plt.axvline(x=maxfreq, linestyle='--')
    plt.show()

A = linfunc(np.ones(2*args.num), spec.freq, peak_freqs)
c, resid, rank, sigma = linalg.lstsq(A, scaled_flux)
lomb_baseline = np.sum(A*c, axis=1) + baseflux.mean()
baseline += lomb_baseline

if args.debug:
    plt.plot(spec.freq, spec.flux,  drawstyle='steps-mid')
    plt.plot(spec.freq, baseline)
    plt.axvline(x=freq0[args.mol], linestyle='--')
#     plt.plot(spec.freq, 0.03*np.sin(np.max(peak_freqs)*spec.freq), 'red')
    plt.show()

spec.flux -= baseline
delv = np.abs(np.average(spec.vel[1:]-spec.vel[:-1]))
n = np.ceil(2*0.4/delv)
rms = np.std(spec.flux[4:-4])*1e3
upper = 3 * np.sqrt(n) * delv * rms
intens = gildas.intens(spec.flux, spec.vel, lim=args.lim)
print('intens = {0[0]} {0[1]} K km/s'.format(intens))
print('rms = {0:.2f} mK, delv = {1:.3f} km/s, upper= {2:.2f} K m/s'.format(rms,
        delv, upper))
mask = [np.abs(spec.vel) <= 20]
plt.plot(spec.vel[mask], spec.flux[mask],  drawstyle='steps-mid')
# plt.axvspan(*args.lim, facecolor='b', alpha=0.5)
plt.axhline(y=0)
plt.axvline(x=0)
plt.savefig('{0}_{1}_baseline.pdf'.format(obsid, args.backend))
plt.show()

np.savetxt(expanduser("~/HssO/Christensen/data/ascii/{}_{:.0f}_{}.dat".format(
                    obsid, freq0[args.mol], args.backend)),
                    np.transpose((spec.vel[mask], spec.flux[mask]*1e3)))
