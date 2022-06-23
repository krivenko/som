# Import HDFArchive and some TRIQS modules
from h5 import HDFArchive
from triqs.gf import *
import triqs.utility.mpi as mpi

# Import main SOM class and utility functions
from som import Som, fill_refreq, compute_tail, reconstruct

n_w = 1001                   # Number of energy slices for the solution
energy_window = (-5.0, 5.0)  # Energy window to search the solution in
tail_max_order = 10          # Maximum tail expansion order to be computed

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
# Right boundary of the histogram, in units of \chi_{min}
acc_params['hist_max'] = 4.0

# Read \chi(i\omega_n) from archive.
# Could be \chi(\tau) or \chi_l as well.
chi_iw = HDFArchive('input.h5', 'r')['chi_iw']

# Set the error bars to a constant (all points of chi_iw are equally important)
error_bars = chi_iw.copy()
error_bars.data[:] = 0.001

# Construct a SOM object
# Norms of spectral functions are known to be 3.0 for both diagonal components
# of chi_iw. It would also be possible to estimate the norms by calling
# som.estimate_boson_corr_spectrum_norms(chi_iw).
cont = Som(chi_iw, error_bars, kind="BosonCorr", norms=[3.0, 3.0])

# Accumulate particular solutions. This may take some time ...
cont.accumulate(**acc_params)

# Construct the final solution as a sum of good particular solutions with
# equal weights. Good particular solutions are those with \chi <= good_d_abs
# and \chi/\chi_{min} <= good_d_rel.
good_chi_abs = 1.0
good_chi_rel = 4.0
cont.compute_final_solution(good_chi_abs=good_chi_abs,
                            good_chi_rel=good_chi_rel,
                            verbosity=1)

# Recover \chi(\omega) on an energy mesh.
# NB: we can use *any* energy window at this point, not necessarily that
# from 'acc_params'.
chi_w = GfReFreq(window=(-5.0, 5.0),
                 n_points=n_w,
                 indices=chi_iw.indices)
fill_refreq(chi_w, cont)

# Do the same, but this time without binning.
chi_w_wo_binning = GfReFreq(window=(-5.0, 5.0),
                            n_points=n_w,
                            indices=chi_iw.indices)
fill_refreq(chi_w_wo_binning, cont, with_binning=False)

# Compute tail coefficients of \chi(\omega)
tail = compute_tail(tail_max_order, cont)

# \chi(i\omega) reconstructed from the solution
chi_rec_iw = chi_iw.copy()
reconstruct(chi_rec_iw, cont)

# On master node, save parameters and results to an archive
if mpi.is_master_node():
    with HDFArchive("results.h5", 'w') as ar:
        ar['acc_params'] = acc_params
        ar['good_chi_abs'] = good_chi_abs
        ar['good_chi_rel'] = good_chi_rel
        ar['chi_iw'] = chi_iw
        ar['chi_rec_iw'] = chi_rec_iw
        ar['chi_w'] = chi_w
        ar['chi_w_wo_binning'] = chi_w_wo_binning
        ar['tail'] = tail
        ar['histograms'] = cont.histograms
