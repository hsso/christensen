#!/usr/bin/python
"""Print Christensen observing log table"""

import numpy as np
from christensen import datadir
from os.path import join

def profile(rhs, rc, m, n, k):
    return ((rhs/rc)**m)*(1+(rhs/rc)**n)**k

skal = 3.9e+28
per = 3.126
dash_par = (2.808, -2.15, 5.093, -4.6142)
solid_par = (5.6, -2.1, 3.2, -3.9)
rh = np.linspace(3.126, 6.5, 100)
dashed = skal*profile(rh, *dash_par)/profile(per, *dash_par)
solid = skal*profile(rh, *solid_par)/profile(per, *solid_par)

np.savetxt(join(datadir, 'ascii', 'dashed.dat'), np.transpose((rh, dashed)))
np.savetxt(join(datadir, 'ascii', 'solid.dat'), np.transpose((rh, solid)))
