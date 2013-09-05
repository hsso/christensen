#!/usr/bin/python

from herschel import HIFISpectrum, hififits
import matplotlib.pyplot as plt
from christensen import datadir
from hsso.gildas import deltadot

obsid = 1342204014

spec = HIFISpectrum(hififits(datadir, obsid, 'WBS', 'H', 'LSB'), 4)
specv = HIFISpectrum(hififits(datadir, obsid, 'WBS', 'V', 'LSB'), 4)
ra = deltadot(spec.start, filename="horizons.txt", column=2).item()
dec = deltadot(spec.start, filename="horizons.txt", column=3).item()
plt.scatter(ra, dec, marker='x')
ra = deltadot(spec.end, filename="horizons.txt", column=2).item()
dec = deltadot(spec.end, filename="horizons.txt", column=3).item()
plt.scatter(ra, dec, marker='x')
plt.scatter(spec.ra, spec.dec)
plt.scatter(specv.ra, specv.dec)
plt.show()
