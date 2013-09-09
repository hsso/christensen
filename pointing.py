#!/usr/bin/python

from herschel import HIFISpectrum, hififits, fwhm
import matplotlib.pyplot as plt
from christensen import datadir
from hsso.gildas import deltadot

obsid = 1342204014

spec = HIFISpectrum(hififits(datadir, obsid, 'WBS', 'H', 'LSB'), 4)
specv = HIFISpectrum(hififits(datadir, obsid, 'WBS', 'V', 'LSB'), 4)
ra_start = deltadot(spec.start, filename="horizons.txt", column=2).item()
dec_start = deltadot(spec.start, filename="horizons.txt", column=3).item()
ra_end = deltadot(spec.end, filename="horizons.txt", column=2).item()
dec_end = deltadot(spec.end, filename="horizons.txt", column=3).item()
ra_mid = deltadot(spec.mid_time, filename="horizons.txt", column=2).item()
dec_mid = deltadot(spec.mid_time, filename="horizons.txt", column=3).item()
plt.plot((ra_start, ra_mid, ra_end), (dec_start, dec_mid, dec_end), marker='x')
beam = fwhm()/3600.
circh = plt.Circle((spec.ra, spec.dec), beam, color='blue', alpha=0.5)
circv = plt.Circle((specv.ra, specv.dec), beam, color='red', alpha=0.5)
ax = plt.gca()
ax.add_patch(circh)
ax.add_patch(circv)
plt.scatter(spec.ra, spec.dec, color='blue')
plt.scatter(specv.ra, specv.dec, color='red')
ax.invert_xaxis()
plt.axis('equal')
plt.show()
