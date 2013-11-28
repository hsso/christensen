#!/usr/bin/python

import numpy as np

ang = 78
alpha = 3*np.pi/2 - ang*np.pi/180
cos, sin = np.cos(alpha), np.sin(alpha)
print("\draw[->,dashed,thick] (axis cs:{0:.3f},{1:.3f}) --\n"
    "(axis cs:{2:.3f},{3:.3f});".format(10*cos, 10*sin, 20*cos, 20*sin))
