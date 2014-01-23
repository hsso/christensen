#!/usr/bin/python

from herschel import HIFISpectrum, hififits, fwhm
import matplotlib.pyplot as plt
from christensen import datadir
from hsso.gildas import deltadot
import numpy as np
import pyfits
from os.path import join

def center(x, y, factor):
    return (x-y)*factor

obsid = 1342204014

hdus = pyfits.open(join(datadir, str(obsid),
    'auxiliary', 'Pointing',
    'hacms1342204014hppv0.1_1373994841254.fits.gz'))

H_s = (260.194933, -47.590554, 268.507438)
A_s = (260.178380, -47.754790, 268.519700)
H_e = (260.189776, -47.590027, 268.524744)
A_e = (260.173150, -47.754260, 268.537060)

spec = HIFISpectrum(hififits(datadir, obsid, 'WBS', 'H', 'LSB'), 4)
specv = HIFISpectrum(hififits(datadir, obsid, 'WBS', 'V', 'LSB'), 4)
ra_start = deltadot(spec.start, filename="horizons.txt", column=2).item()
dec_start = deltadot(spec.start, filename="horizons.txt", column=3).item()
ra_end = deltadot(spec.end, filename="horizons.txt", column=2).item()
dec_end = deltadot(spec.end, filename="horizons.txt", column=3).item()
ra_mid = deltadot(spec.mid_time, filename="horizons.txt", column=2).item()
dec_mid = deltadot(spec.mid_time, filename="horizons.txt", column=3).item()
# RA and Dec of synthetic beam
ra_beam = np.average((spec.ra, specv.ra))
dec_beam = np.average((spec.dec, specv.dec))
spec.ra -= ra_beam
specv.ra -= ra_beam
spec.dec -= dec_beam
specv.dec -= dec_beam
spec.ra *= 3600
spec.dec *= 3600
specv.ra *= 3600
specv.dec *= 3600
ra_start = center(ra_start, ra_beam, 3600)
dec_start = center(dec_start, dec_beam, 3600)
ra_end = center(ra_end, ra_beam, 3600)
dec_end = center(dec_end, dec_beam, 3600)
ra_mid = center(ra_mid, ra_beam, 3600)
dec_mid = center(dec_mid, dec_beam, 3600)

print(ra_beam, dec_beam)
h_start = center(np.array(H_s[:-1]), np.array((ra_beam, dec_beam)), 3600)
h_end = center(np.array(H_e[:-1]), np.array((ra_beam, dec_beam)), 3600)
a_start = center(np.array(A_s[:-1]), np.array((ra_beam, dec_beam)), 3600)
a_end = center(np.array(A_e[:-1]), np.array((ra_beam, dec_beam)), 3600)

print(spec.ra, spec.dec)
print(specv.ra, specv.dec)
print(ra_start, dec_start)
print(ra_mid, dec_mid)
print(ra_end, dec_end)
beam = fwhm()/2.
print(beam)
print(np.linalg.norm(np.array((ra_mid, dec_mid))-
        np.array((spec.ra, spec.dec))))
print(np.linalg.norm(np.array((ra_mid, dec_mid))-
        np.array((specv.ra, specv.dec))))
circh = plt.Circle((spec.ra, spec.dec), beam, color='blue', alpha=0.5)
circv = plt.Circle((specv.ra, specv.dec), beam, color='red', alpha=0.5)
ax = plt.gca()
ax.arrow(ra_start, dec_start, ra_end-ra_start, dec_end-dec_start,
    fc='k', ec='k', head_width=1, head_length=1)
ax.add_patch(circh)
ax.add_patch(circv)
plt.scatter(spec.ra, spec.dec, color='blue')
plt.scatter(specv.ra, specv.dec, color='red')
ax.arrow(h_start[0], h_start[1], *(h_end-h_start),
    head_width=1, head_length=1, color='green', ec='g')
ax.invert_xaxis()
plt.axis('equal')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.savefig('pointing.pdf')
plt.show()
