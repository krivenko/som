# Import HDFArchive and some TRIQS modules
from h5 import HDFArchive
from triqs.gf import *
import triqs.utility.mpi as mpi

# Import main SOM class and utility functions
from som import (Som,
                 fill_refreq,
                 compute_tail,
                 reconstruct,
                 estimate_boson_corr_spectrum_norms)

n_w = 501                    # Number of energy slices for the solution
energy_window = (-4.0, 4.0)  # Energy window to search the solution in
tail_max_order = 10          # Maximum tail expansion order to be computed

# Parameters for Som.accumulate()
acc_params = {'energy_window': energy_window}
acc_params['verbosity'] = 2  # Verbosity level
acc_params['l'] = 2000       # Number of particular solutions to accumulate
acc_params['f'] = 100        # Number of global updates
acc_params['t'] = 50         # Number of local updates per global update

# Read \chi(i\Omega_n), \chi(\tau) and \chi_l from an archive.
with HDFArchive('input.h5', 'r') as ar:
    chi_iw = ar['chi_iw']
    chi_tau = ar['chi_tau']
    chi_l = ar['chi_l']

#
# Analytically continue \chi(i\Omega_n), \chi(\tau) and \chi_l
#

for mesh_name, chi in (("ImFreq", chi_iw),
                       ("ImTime", chi_tau),
                       ("Legendre", chi_l)):

    # Estimate norms of spectral functions from input data
    norms = estimate_boson_corr_spectrum_norms(chi)

    # Set the error bars to a constant
    error_bars = chi.copy()
    error_bars.data[:] = 0.001

    # Construct a SOM object
    cont = Som(chi, error_bars, kind="BosonAutoCorr", norms=norms)

    # Accumulate particular solutions
    cont.accumulate(**acc_params)

    # Construct the final solution
    cont.compute_final_solution(good_chi_rel=4.0, verbosity=1)

    # Recover \chi(\omega) on an energy mesh.
    chi_w = GfReFreq(window=energy_window, n_points=n_w, indices=chi.indices)
    fill_refreq(chi_w, cont)

    # \chi reconstructed from the solution
    chi_rec = chi.copy()
    reconstruct(chi_rec, cont)

    # On master node, save results to an archive
    if mpi.is_master_node():
        with HDFArchive("results.h5", 'a') as ar:
            ar.create_group(mesh_name)
            gr = ar[mesh_name]
            gr['norms'] = norms
            gr['chi'] = chi
            gr['chi_w'] = chi_w
            gr['chi_rec'] = chi_rec
