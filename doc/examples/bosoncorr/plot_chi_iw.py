from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

# Plot input and reconstructed \chi(i\omega_n), real parts
oplot(ar['chi_iw'][0,0],     mode='R', linewidth=0.8, label="$\chi'_0(i\\omega_n)$")
oplot(ar['chi_rec_iw'][0,0], mode='R', linewidth=0.8, label="$\chi'_\mathrm{0,rec}(i\\omega_n)$")
oplot(ar['chi_iw'][1,1],     mode='R', linewidth=0.8, label="$\chi'_1(i\\omega_n)$")
oplot(ar['chi_rec_iw'][1,1], mode='R', linewidth=0.8, label="$\chi'_\mathrm{1,rec}(i\\omega_n)$")

plt.xlim((0, 20))
plt.ylabel("$\chi'(i\\omega_n)$")
plt.legend(loc="upper center")
plt.show()

# Plot input and reconstructed \chi(i\omega_n), imaginary parts
oplot(ar['chi_iw'][0,0],     mode='I', linewidth=0.8, label="$\chi''_0(i\\omega_n)$")
oplot(ar['chi_rec_iw'][0,0], mode='I', linewidth=0.8, label="$\chi''_\mathrm{0,rec}(i\\omega_n)$")
oplot(ar['chi_iw'][1,1],     mode='I', linewidth=0.8, label="$\chi''_1(i\\omega_n)$")
oplot(ar['chi_rec_iw'][1,1], mode='I', linewidth=0.8, label="$\chi''_\mathrm{1,rec}(i\\omega_n)$")

plt.xlim((0, 20))
plt.ylabel("$\chi(i\\omega_n)$")
plt.legend(loc="lower center")
