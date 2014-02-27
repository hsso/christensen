#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser, join
import argparse
from scipy import linalg, fftpack
from hsso import gildas
from hsso.class_utils import pgram_peaks, linfunc, fitfunc
from herschel import HIFISpectrum, hififits
from christensen import datadir, freq0

def fileout(suffix):
    """Return full path to data file"""
    filename = "{}_{}{}.dat".format(obsid, args.backend, suffix)
    return join(datadir, 'ascii', filename)

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--backend', default='WBS', choices=('HRS', 'WBS'))
parser.add_argument('--sideband', default='', choices=('', 'LSB', 'USB'))
parser.add_argument('--subband', default=0, type=int, choices=range(5))
parser.add_argument('-d', '--debug', action='store_true', help='debug mode')
parser.add_argument('--fftlim', default=0, type=float,
                    help='FFT high frequency limit')
parser.add_argument('-n', '--num', default=8, type=int,
                    help='number of sine waves')
parser.add_argument('-m', '--mol', default='H2O', choices=('H2O', 'NH3'))
parser.add_argument('--deg', default=2, type=int,
                    help='baseline polynomial degree')
parser.add_argument('--twiny', action='store_true')
parser.add_argument('--fold', action='store_true')
parser.add_argument("--lim", nargs=2, type=float, default=(-1, 1),
                    help='line limits in km/s')
args = parser.parse_args()

subband = {'HRS': 1, 'WBS': 4}
sideband = {'H2O': 'LSB', 'NH3': 'USB'}
fftlim = {'HRS': 1e3, 'WBS': 2.5e2}
if not args.subband: args.subband = subband[args.backend]
if not args.sideband: args.sideband = sideband[args.mol]
if not args.fftlim: args.fftlim = fftlim[args.backend]

obsid = 1342204014
spec = HIFISpectrum(hififits(datadir, obsid, args.backend, 'H', args.sideband),
                    args.subband, freq0=freq0[args.mol])
specv = HIFISpectrum(hififits(datadir, obsid, args.backend, 'V', args.sideband),
                    args.subband, freq0=freq0[args.mol])

if args.mol == "NH3":
    spec.fold()
    specv.fold()
    spec.add(specv)
    spec.fftbase(args.fftlim, line=(), plot=args.debug)
    if args.backend == "HRS": spec.resample()
    spec.save(fileout('_'+args.mol), "fluxcal")
    if args.debug: spec.plot()
elif args.fold:
    spec.fold()
    spec.fftbase(args.fftlim, shift=-0., linelim=1., plot=args.debug)
    if args.backend == "HRS": spec.resample()
    spec.save(fileout('-H_fluxcal'), "fluxcal")
    if args.debug: spec.plot()
    print(spec.intens, spec.error, spec.snr)

    specv.fold()
    specv.fftbase(args.fftlim, plot=args.debug)
    if args.backend == "HRS": specv.resample()
    specv.save(fileout('-V_fluxcal'), "fluxcal")
    if args.debug: specv.plot()
    print(specv.intens, specv.error, specv.snr)

    spec.add(specv)
    spec.fftbase(args.fftlim, shift=-0., plot=args.debug)
    spec.save(fileout('_fluxcal'), "fluxcal")
    if args.debug: spec.plot()
    print(spec.intens, spec.error, spec.snr)
else:
    spec.scale((-60, 10))
    spec.save(fileout('-H'))
    spec.fftbase(args.fftlim, shift=-0.4, linelim=1., throw=True, plot=args.debug)
    spec.save(fileout('-H_unfolded'), "baseline")
    spec.save(fileout('-H_unfold_fluxcal'), "fluxcal")
    if args.debug: spec.plot()

    specv.scale((-60, 10))
    specv.save(fileout('-V'))
    specv.fftbase(args.fftlim, linelim=1., throw=True, plot=args.debug)
    specv.save(fileout('-V_unfolded'), "baseline")
    specv.save(fileout('-V_unfold_fluxcal'), "fluxcal")
    if args.debug: specv.plot()

    spec.add(specv)
    spec.save(fileout('_aver'))
    spec.fftbase(args.fftlim, shift=-0.2, linelim=1., throw=True, plot=args.debug)
    spec.save(fileout('_unfolded'), "baseline")
    spec.save(fileout('_unfold_fluxcal'), "fluxcal")
    if args.debug: spec.plot()

# # Lomb-Scargle periodogram
# spec.baseflux -= spec.baseline
# scaled_flux = spec.baseflux-spec.baseflux.mean()
# # frequency
# f = np.linspace(1e2, 2e4, 1e4)
# pgram, peak_freqs, peak_flux = pgram_peaks(spec.freq, scaled_flux, f, args.num)
# if args.debug:
#     plt.loglog(f, pgram)
#     for maxfreq in peak_freqs:
#         plt.axvline(x=maxfreq, linestyle='--')
#     plt.show()
# 
# A = linfunc(np.ones(2*args.num), spec.freq, peak_freqs)
# c, resid, rank, sigma = linalg.lstsq(A, scaled_flux)
# lomb_baseline = np.sum(A*c, axis=1) + spec.baseflux.mean()
# spec.baseline += lomb_baseline
# 
# if args.debug:
#     plt.plot(spec.freq, spec.flux,  drawstyle='steps-mid')
#     plt.plot(spec.freq, spec.baseline)
#     plt.axvline(x=freq0[args.mol], linestyle='--')
# #     plt.plot(spec.freq, 0.03*np.sin(np.max(peak_freqs)*spec.freq), 'red')
#     plt.show()
# 
# spec.flux -= spec.baseline
# delv = np.abs(np.average(spec.vel[1:]-spec.vel[:-1]))
# n = np.ceil(2*0.4/delv)
# rms = np.std(spec.flux[4:-4])*1e3
# upper = 3 * np.sqrt(n) * delv * rms
# intens = gildas.intens(spec.flux, spec.vel, lim=args.lim)
# print('intens = {0[0]} {0[1]} K km/s'.format(intens))
# print('rms = {0:.2f} mK, delv = {1:.3f} km/s, upper= {2:.2f} K m/s'.format(rms,
#         delv, upper))
# mask = [np.abs(spec.vel) <= 20]
# plt.plot(spec.vel[mask], spec.flux[mask],  drawstyle='steps-mid')
# # plt.axvspan(*args.lim, facecolor='b', alpha=0.5)
# plt.axhline(y=0)
# plt.axvline(x=0)
# plt.savefig('{0}_{1}_baseline.pdf'.format(obsid, args.backend))
# plt.show()
# 
# np.savetxt(expanduser("~/HssO/Christensen/data/ascii/{}_{:.0f}_{}.dat".format(
#                     obsid, freq0[args.mol], args.backend)),
#                     np.transpose((spec.vel[mask], spec.flux[mask])))
