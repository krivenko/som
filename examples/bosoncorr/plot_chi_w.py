from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from matplotlib import pyplot as plt
from pytriqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

# Plot imaginary part of the susceptibility on the real axis
oplot(ar['chi_w'][0,0], mode='I', linewidth=0.8, label="$\chi''_0(\\omega)$")
oplot(ar['chi_w'][1,1], mode='I', linewidth=0.8, label="$\chi''_1(\\omega)$")

plt.xlim((-5.0,5.0))
plt.ylim((-1.5,1.5))
plt.ylabel("$\chi(\\omega)$")
plt.legend(loc = "lower right")
