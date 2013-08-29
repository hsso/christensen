#!/usr/bin/python

from christensen import obsids, datadir, paperdir
from os.path import join
import glob
import pyfits
from datetime import datetime
from hsso.gildas import frac_day, deltadot

def hierarch_key(header, value):
    for i in header.ascardlist().keys():
        if i.find('META_') > 0 and header[i] == value:
            return header[i[4:]]

f = open(join(paperdir, 'log.txt'), 'w')
for obsid in obsids:
    hdus = pyfits.open(glob.glob(join(datadir, str(obsid), '*fits.gz'))[0])
    date_obs = hdus[0].header['DATE-OBS']
    date_end = hdus[0].header['DATE-END']
    start = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S.%f")
    end = datetime.strptime(date_end, "%Y-%m-%dT%H:%M:%S.%f")
    exp = end-start
    mid_time = start + exp/2
    exp_min = exp.seconds/60.
    ins = hdus[0].header['INSTRUME']
    if ins == "PACS":
        scan_ang = "{:.0f}".format(hierarch_key(hdus[0].header, 'mapScanAngle'))
    else:
        scan_ang = ""
    # calculate phase angle
    rh = deltadot(mid_time, filename="horizons.txt", column=4).item()
    delta = deltadot(mid_time, filename="horizons.txt", column=6).item()
    phi = deltadot(mid_time, filename="horizons.txt", column=8).item()
    f.write("{:%Y-%m}-{:06.3f} & {} & {} & {:.1f} & {} & {:.2f} & {:.2f} & {:.2f}\\\\\n".format(mid_time,
        frac_day(mid_time), ins, obsid, exp_min, scan_ang, rh, delta, phi))
