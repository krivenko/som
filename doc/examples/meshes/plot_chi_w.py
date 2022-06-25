from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

#
# Plot susceptibility as a function of real frequencies
#

# Real parts
oplot(ar['ImFreq']['chi_w'][0, 0],
      mode='R', lw=0.8, label=r"$\chi'(\omega)$ from $\chi(i\Omega_n)$")
oplot(ar['ImTime']['chi_w'][0, 0],
      mode='R', lw=0.8, label=r"$\chi'(\omega)$ from $\chi(\tau)$")
oplot(ar['Legendre']['chi_w'][0, 0],
      mode='R', lw=0.8, label=r"$\chi'(\omega)$ from $\chi(\ell)$")

# Imaginary parts
oplot(ar['ImFreq']['chi_w'][0, 0],
      mode='I', lw=0.8, label=r"$\chi''(\omega)$ from $\chi(i\Omega_n)$")
oplot(ar['ImTime']['chi_w'][0, 0],
      mode='I', lw=0.8, label=r"$\chi''(\omega)$ from $\chi(\tau)$")
oplot(ar['Legendre']['chi_w'][0, 0],
      mode='I', lw=0.8, label=r"$\chi''(\omega)$ from $\chi(\ell)$")

# Spectrum normalization constants extracted from input data
norms_text = [
      r"Spectrum norm from $\chi(i\Omega_n)$: %.3f" % ar['ImFreq']['norms'][0],
      r"Spectrum norm from $\chi(\tau)$: %.3f" % ar['ImTime']['norms'][0],
      r"Spectrum norm from $\chi(\ell)$: %.3f" % ar['Legendre']['norms'][0]
]
norms_text = "\n".join(norms_text)
plt.text(0, -1.5, norms_text)

plt.ylabel(r"$\chi(\omega)$")
plt.legend(loc="upper left")
