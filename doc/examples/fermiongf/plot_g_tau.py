from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

g_tau = ar['g_tau']
g_rec_tau = ar['g_rec_tau']

# Plot input and reconstructed G(\tau)
oplot(g_tau[0, 0],     mode='R', lw=0.8, label=r"$G_{00}(\tau)$")
oplot(g_tau[1, 1],     mode='R', lw=0.8, label=r"$G_{11}(\tau)$")
oplot(g_rec_tau[0, 0], mode='R', lw=0.8, label=r"$G_{00}^\mathrm{rec}(\tau)$")
oplot(g_rec_tau[1, 1], mode='R', lw=0.8, label=r"$G_{11}^\mathrm{rec}(\tau)$")

plt.xlim((0, g_tau.mesh.beta))
plt.ylabel(r"$G(\tau)$")
plt.legend(loc="lower center")
