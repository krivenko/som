from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from matplotlib import pyplot as plt
from pytriqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('example.h5', 'r')

# Plot the spectral function
oplot(ar['g_w'], mode='S', label = "$A(\\epsilon)$")

plt.show()
