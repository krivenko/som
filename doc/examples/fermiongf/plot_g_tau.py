from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

# Plot input and reconstructed G(\tau)
oplot(ar['g_tau'][0,0],     mode='R', linewidth=0.8, label="$G_{00}(\\tau)$")
oplot(ar['g_tau'][1,1],     mode='R', linewidth=0.8, label="$G_{11}(\\tau)$")
oplot(ar['g_rec_tau'][0,0], mode='R', linewidth=0.8, label="$G_{00}^\mathrm{rec}(\\tau)$")
oplot(ar['g_rec_tau'][1,1], mode='R', linewidth=0.8, label="$G_{11}^\mathrm{rec}(\\tau)$")

plt.xlim((0, 20))
plt.ylabel("$G(\\tau)$")
plt.legend(loc="lower center")
