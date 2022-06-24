from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot

# Read data from archive
ar = HDFArchive('results.h5', 'r')

# Plot spectral functions
oplot(ar['error_bars']['g_w'][0, 0],
      mode='S', lw=0.8, label=r"Error bars")
oplot(ar['cov_matrix']['g_w'][0, 0],
      mode='S', lw=0.8, label=r"Covariance matrix")
oplot(ar['cov_matrix_fl']['g_w'][0, 0],
      mode='S', lw=0.8, label=r"Covariance matrix with filtering")

plt.ylim((0, 0.4))
plt.ylabel(r"$A(\omega)$")
plt.legend(loc="upper left")
