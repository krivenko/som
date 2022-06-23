from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

chi_w = ar['chi_w']
chi_w_wob = ar['chi_w_wo_binning']
tail = ar['tail']

# Plot imaginary part of the correlator on the real axis
# with and without binning
oplot(chi_w_wob[0, 0],
      mode='I', lw=0.8, label=r"$\chi''_0(\omega)$, w/o binning")
oplot(chi_w_wob[1, 1],
      mode='I', lw=0.8, label=r"$\chi''_1(\omega)$, w/o binning")
oplot(chi_w[0, 0], mode='I', lw=0.8, label=r"$\chi''_0(\omega)$")
oplot(chi_w[1, 1], mode='I', lw=0.8, label=r"$\chi''_1(\omega)$")

plt.xlim((-2.5, 2.5))
plt.ylim((-0.15, 0.15))
plt.ylabel(r"$\chi''(\omega)$")
plt.legend(loc="lower right")
