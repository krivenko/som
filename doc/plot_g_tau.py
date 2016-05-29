from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from matplotlib import pyplot as plt
from pytriqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('example.h5', 'r')

# Plot input and reconstructed G(\tau)
oplot(ar['g_tau'],     mode='R', label = "$G(\\tau)$")
oplot(ar['g_rec_tau'], mode='R', label = "$G_\mathrm{rec}(\\tau)$")

ax = plt.gca()
ax.set_ylabel('')
ax.legend(loc = 'lower center')

plt.show()
