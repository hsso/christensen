#!/usr/bin/python

import tarfile
from os.path import join
from christensen import datadir, hifi_obsid

with tarfile.open("christensen_updp.tar.gz", "w:gz") as tar:
    for backend in ['WBS', 'HRS']:
        outfile = "{}_{}.pdf".format(hifi_obsid, backend)
        source_file = join(datadir, outfile)
        tar.add(source_file, arcname=join('Christensen', 'HFI', outfile))
tar.close()
