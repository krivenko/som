from h5 import HDFArchive
from triqs.stat import Histogram
from matplotlib import pyplot as plt
from som import count_good_solutions
import numpy

# Read data from archive
ar = HDFArchive('results.h5', 'r')

hist = ar['histograms'][0]
good_chi_abs = ar['good_chi_abs']
good_chi_rel = ar['good_chi_rel']

# Plot histogram
chi_grid = numpy.array([hist.mesh_point(n) for n in range(len(hist))])
plt.plot(chi_grid, hist.data, label=r"Histogram for $A(\omega)$")

plt.xlim((chi_grid[0], chi_grid[-1]))
plt.xlabel(r"$\chi$")
plt.ylabel(r"$P(\chi)$")
plt.legend(loc="upper right")

# Count good solutions using saved histograms
n_good_sol = count_good_solutions(hist,
                                  good_chi_abs=good_chi_abs,
                                  good_chi_rel=good_chi_rel)

plt.text(0.3125, 90, "# good solutions in $A(\\omega)$: %d" % n_good_sol)
