from pytriqs.gf import *
from pytriqs.archive import HDFArchive
from matplotlib import pyplot as plt
from pytriqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

# Plot input and reconstructed \chi(i\omega_n)
oplot(ar['chi_iw'][0,0],     mode='R', linewidth=0.8, label="$\chi_0(i\\omega_n)$")
oplot(ar['chi_rec_iw'][0,0], mode='R', linewidth=0.8, label="$\chi_\mathrm{0,rec}(i\\omega_n)$")
oplot(ar['chi_iw'][1,1],     mode='R', linewidth=0.8, label="$\chi_1(i\\omega_n)$")
oplot(ar['chi_rec_iw'][1,1], mode='R', linewidth=0.8, label="$\chi_\mathrm{1,rec}(i\\omega_n)$")

plt.xlim((0, 3))
plt.ylabel("$\chi'(i\\omega_n)$")
plt.legend(loc="upper right")
