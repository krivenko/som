from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

# Plot input and reconstructed \chi(i\Omega_n)
oplot(ar['ImFreq']['chi'][0, 0], mode='R', lw=0.8, label=r"$\chi(i\Omega_n)$")
oplot(ar['ImFreq']['chi_rec'][0, 0], mode='R', lw=0.8,
      label=r"$\chi^{rec}(i\Omega_n)$")

plt.ylabel(r"$\chi(i\Omega_n)$")
plt.legend(loc="upper right")
plt.show()

# Plot input and reconstructed \chi(\tau)
oplot(ar['ImTime']['chi'][0, 0], mode='R', lw=0.8, label=r"$\chi(\tau)$")
oplot(ar['ImTime']['chi_rec'][0, 0], mode='R', lw=0.8,
      label=r"$\chi^{rec}(\tau)$")

plt.ylabel(r"$\chi(\tau)$")
plt.legend(loc="upper center")
plt.show()

# Plot input and reconstructed \chi(\ell)
oplot(ar['Legendre']['chi'][0, 0], mode='R', lw=0.8, label=r"$\chi(\ell)$")
oplot(ar['Legendre']['chi_rec'][0, 0], mode='R', lw=0.8,
      label=r"$\chi^{rec}(\ell)$")

plt.ylabel(r"$\chi(\ell)$")
plt.legend(loc="upper right")
plt.show()
