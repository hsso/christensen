#!/usr/bin/python

import tarfile
from os.path import join
from christensen import datadir, figsdir, hifi_obsid, obsids

with tarfile.open("christensen_updp.tar.gz", "w:gz") as tar:
    # HIFI data
    for backend in ['WBS', 'HRS']:
        for ext in ['pdf', 'fits']:
            outfile = "{}_{}.{}".format(hifi_obsid, backend, ext)
            source_file = join(datadir, outfile)
            tar.add(source_file, arcname=join('Christensen', 'HIFI', outfile))
    # PACS data
    for obsid in obsids[::2]: # First obsid in each pair
        for band in ['blue', 'red']:
            outfile = "{}_{}.png".format(obsid, band)
            tar.add(source_file, arcname=join('Christensen', 'PACS', outfile))

tar.close()
