from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

chi_iw = ar['chi_iw']
chi_rec_iw = ar['chi_rec_iw']

# Plot real part of input and reconstructed \chi(i\omega_n)
oplot(chi_iw[0, 0],     mode='R', lw=0.8, label=r"$\chi'_0(i\omega_n)$")
oplot(chi_rec_iw[0, 0], mode='R', lw=0.8, label=r"$\chi'_{0,rec}(i\omega_n)$")
oplot(chi_iw[1, 1],     mode='R', lw=0.8, label=r"$\chi'_1(i\omega_n)$")
oplot(chi_rec_iw[1, 1], mode='R', lw=0.8, label=r"$\chi'_{1,rec}(i\omega_n)$")

plt.xlim((0, 20.0))
plt.ylabel(r"$\chi'(i\omega_n)$")
plt.legend(loc="upper center")
plt.show()

# Plot imaginary part of input and reconstructed \chi(i\omega_n)
oplot(chi_iw[0, 0],     mode='I', lw=0.8, label=r"$\chi''_0(i\omega_n)$")
oplot(chi_rec_iw[0, 0], mode='I', lw=0.8, label=r"$\chi''_{0,rec}(i\omega_n)$")
oplot(chi_iw[1, 1],     mode='I', lw=0.8, label=r"$\chi''_1(i\omega_n)$")
oplot(chi_rec_iw[1, 1], mode='I', lw=0.8, label=r"$\chi''_{1,rec}(i\omega_n)$")

plt.xlim((0, 20.0))
plt.ylabel(r"$\chi''(i\omega_n)$")
plt.legend(loc="lower center")
