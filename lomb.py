#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg

def pgram_peaks(freq, flux, f, num):
    from scipy.signal import lombscargle
    normval = freq.shape[0]
    pgram = lombscargle(freq, flux, f)/normval
    pgram = np.sqrt(4*(pgram/normval))
    # find local maxima not at the edge of the periodogram
    maxmask = np.r_[False, pgram[1:] > pgram[:-1]] &\
                np.r_[pgram[:-1] > pgram[1:], False]
    sortarg = np.argsort(pgram[maxmask])
    peak_freqs = f[maxmask][sortarg[-num:]]
    peak_flux = pgram[maxmask][sortarg[-num:]]
    return pgram, peak_freqs, peak_flux

def linfunc(a, x, peak_freqs):
    """Target function"""
    num = len(peak_freqs)
    sinwave = [a[i]*np.sin(peak_freqs[i]*x) for i in range(len(peak_freqs))]
    coswave = [a[num+i]*np.cos(peak_freqs[i]*x) for i in range(len(peak_freqs))]
    return np.column_stack(sinwave+coswave)

freq, flux = np.loadtxt('1342204014.txt', usecols=(0,2), unpack=True)
flux -= flux.mean()
num = 40

# Lomb-Scargle periodogram
# frequency
f = np.linspace(1e1, 3e4, 1e4)
pgram, peak_freqs, peak_flux = pgram_peaks(freq, flux,
        f, num)
plt.loglog(f, pgram)
for maxfreq in peak_freqs:
    plt.axvline(x=maxfreq, linestyle='--')
plt.show()

A = linfunc(np.ones(2*num), freq, peak_freqs)
c, resid, rank, sigma = linalg.lstsq(A, flux)
lomb_baseline = np.sum(A*c, axis=1)
plt.plot(freq, flux)
plt.plot(freq, lomb_baseline)
plt.show()

