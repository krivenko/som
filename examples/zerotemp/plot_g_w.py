from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from matplotlib import pyplot as plt
from pytriqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

# Plot the spectral function
oplot(ar['g_w'], mode='S', linewidth=0.8, label="$A(\\omega)$")

plt.xlim((0,10.0))
plt.ylim((0,0.3))
plt.ylabel("$A(\\omega)$")
plt.legend(loc = "upper right")
