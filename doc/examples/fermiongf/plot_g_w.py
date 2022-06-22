from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from triqs.plot.mpl_interface import oplot
import numpy

# Read data from archive
ar = HDFArchive('results.h5', 'r')

g_w = ar['g_w']
g_w_wob = ar['g_w_wo_binning']
tail = ar['tail']

# Plot spectral functions with and without binning
oplot(g_w_wob[0, 0], mode='S', lw=0.8, label=r"$A_0(\omega)$, w/o binning")
oplot(g_w_wob[1, 1], mode='S', lw=0.8, label=r"$A_1(\omega)$, w/o binning")
oplot(g_w[0, 0], mode='S', linewidth=0.8, label=r"$A_0(\omega)$")
oplot(g_w[1, 1], mode='S', linewidth=0.8, label=r"$A_1(\omega)$")

plt.ylim((0, 0.4))
plt.ylabel(r"$A(\omega)$")
plt.legend(loc="upper left")

# Plot tail coefficients
tail_orders = numpy.array(range(tail.shape[0]))
tail_ax = inset_axes(plt.gca(), width="30%", height="60%", loc="upper right")
tail_ax.semilogy(tail_orders, numpy.abs(tail[:, 0, 0]), label=r"$G_0(\omega)$")
tail_ax.semilogy(tail_orders, numpy.abs(tail[:, 1, 1]), label=r"$G_1(\omega)$")
tail_ax.set_xlim((tail_orders[0], tail_orders[-1]))
tail_ax.set_ylabel("$|a_i|$")
tail_ax.yaxis.tick_right()
tail_ax.legend()
