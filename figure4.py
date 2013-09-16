#!/usr/bin/python
"""Print Christensen observing log table"""

rc=2.808
m=-2.15
n=5.093
k=-4.6142
skal=3.9e+28
per=3.126
alfg=skal/(((per/rc)**m)*(1+(per/rc)**n)**k)
g=alfg*(((rhs/rc)**m)*(1+(rhs/rc)**n)**k)

# solid curve:
rcs=5.6
ms=-2.1
ns=3.2
ks=-3.9
skal=3.9e+28
alfs=skal/(((per/rcs)**ms)*(1+(per/rcs)**ns)**ks)
gs=alfs*(((rhs/rcs)**ms)*(1+(rhs/rcs)**ns)**ks)
