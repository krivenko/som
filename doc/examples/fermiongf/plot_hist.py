from h5 import HDFArchive
from triqs.stat import Histogram
from matplotlib import pyplot as plt
from som import count_good_solutions
import numpy

# Read data from archive
ar = HDFArchive('results.h5', 'r')

hists = ar['histograms']
good_chi_abs = ar['good_chi_abs']
good_chi_rel = ar['good_chi_rel']

# Plot histograms
chi_min, chi_max = numpy.inf, 0
for n, hist in enumerate(hists):
    chi_grid = numpy.array([hist.mesh_point(n) for n in range(len(hist))])
    chi_min = min(chi_min, chi_grid[0])
    chi_max = max(chi_min, chi_grid[-1])
    plt.plot(chi_grid, hist.data, label=r"Histogram for $A_%d(\omega)$" % n)

plt.xlim((chi_min, chi_max))
plt.xlabel(r"$\chi$")
plt.ylabel(r"$P(\chi)$")
plt.legend(loc="upper right")

# Count good solutions using saved histograms
n_good_sol = [count_good_solutions(h,
                                   good_chi_abs=good_chi_abs,
                                   good_chi_rel=good_chi_rel)
              for h in hists]

plt.text(0.06, 33,
         ("# good solutions in $A_0(\\omega)$: %d\n"
          "# good solutions in $A_1(\\omega)$: %d") %
         (n_good_sol[0], n_good_sol[1]))
