from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from matplotlib import pyplot as plt
from pytriqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

# Plot the spectral functions
oplot(ar['g_w'][0,0], mode='S', linewidth=0.8, label="$A_0(\\omega)$")
oplot(ar['g_w'][1,1], mode='S', linewidth=0.8, label="$A_1(\\omega)$")

plt.ylim((0,0.4))
plt.ylabel("$A(\\omega)$")
