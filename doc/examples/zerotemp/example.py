# Import HDFArchive and some TRIQS modules
from h5 import HDFArchive
from triqs.gf import *
import triqs.utility.mpi as mpi

# Import main SOM class and utility functions
from som import Som, fill_refreq, compute_tail, reconstruct

n_w = 1001                 # Number of energy slices for the solution
energy_window = (0, 10.0)  # Energy window to search the solution in
tail_max_order = 10        # Maximum tail expansion order to be computed

# Parameters for Som.accumulate()
acc_params = {'energy_window': energy_window}
# Verbosity level
acc_params['verbosity'] = 2
# Number of particular solutions to accumulate
acc_params['l'] = 2000
# Number of global updates
acc_params['f'] = 100
# Number of local updates per global update
acc_params['t'] = 50
# Accumulate histogram of the objective function values
acc_params['make_histograms'] = True

# Read G(\tau) from archive.
# Could be G(i\omega_n) or G_l as well.
g_tau = HDFArchive('example.h5', 'r')['g_tau']

# Set the error bars to a constant (all points of g_tau are equally important)
error_bars = g_tau.copy()
error_bars.data[:] = 0.01

# Construct a SOM object
# Norm of the spectral function is known to be 0.5.
cont = Som(g_tau, error_bars, kind="ZeroTemp", norms=[0.5])

# Accumulate particular solutions. This may take some time ...
cont.accumulate(**acc_params)

# Construct the final solution as a sum of good particular solutions with
# equal weights. Good particular solutions are those with \chi <= good_d_abs
# and \chi/\chi_{min} <= good_d_rel.
good_chi_abs = 0.002
good_chi_rel = 2.0
cont.compute_final_solution(good_chi_abs=good_chi_abs,
                            good_chi_rel=good_chi_rel,
                            verbosity=1)

# Recover G(\omega) on an energy mesh.
# NB: we can use *any* energy window at this point, not necessarily that
# from 'acc_params'.
g_w = GfReFreq(window=(0, 10.0),
               n_points=n_w,
               indices=[0])
fill_refreq(g_w, cont)

# Do the same, but this time without binning.
g_w_wo_binning = GfReFreq(window=(0, 10.0),
                          n_points=n_w,
                          indices=g_tau.indices)
fill_refreq(g_w_wo_binning, cont, with_binning=False)

# Compute tail coefficients of G(\omega)
tail = compute_tail(tail_max_order, cont)

# G(\tau) reconstructed from the solution
g_rec_tau = g_tau.copy()
reconstruct(g_rec_tau, cont)

# On master node, save parameters and results to an archive
if mpi.is_master_node():
    with HDFArchive("results.h5", 'w') as ar:
        ar['acc_params'] = acc_params
        ar['good_chi_abs'] = good_chi_abs
        ar['good_chi_rel'] = good_chi_rel
        ar['g_tau'] = g_tau
        ar['g_rec_tau'] = g_rec_tau
        ar['g_w'] = g_w
        ar['g_w_wo_binning'] = g_w_wo_binning
        ar['tail'] = tail
        ar['histograms'] = cont.histograms
