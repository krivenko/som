# Import HDFArchive and some TRIQS modules
from h5 import HDFArchive
from triqs.gf import *
import triqs.utility.mpi as mpi

# Import main SOM class
from som import Som
# Import statistical analysis functions
from som.spectral_stats import spectral_avg, spectral_disp, spectral_corr

energy_window = (-4.0, 4.0)  # Energy window to search the solution in

# Parameters for Som.accumulate()
acc_params = {'energy_window': energy_window}
acc_params['verbosity'] = 2  # Verbosity level
acc_params['l'] = 1000       # Number of particular solutions to accumulate
acc_params['f'] = 100        # Number of global updates
acc_params['t'] = 50         # Number of local updates per global update

# Read G(\tau) from an archive.
g_tau = HDFArchive('input.h5', 'r')['g_tau']

# Set the error bars to a constant (all points of g_tau are equally important)
error_bars = g_tau.copy()
error_bars.data[:] = 0.01

# Construct a SOM object
cont = Som(g_tau, error_bars, kind="FermionGf")

# Accumulate particular solutions. This may take some time ...
cont.accumulate(**acc_params)

#
# Compute statistical characteristics of the ensemble of the accumulated
# particular solutions.
#

# We use a set of energy intervals centered around points of this mesh
w_mesh = MeshReFreq(*energy_window, 401)

for i in (0, 1):
    # Compute spectral averages using different resolution functions.
    avg_rect = spectral_avg(cont, i, w_mesh, "rectangle")
    avg_lorentz = spectral_avg(cont, i, w_mesh, "lorentzian")
    avg_gauss = spectral_avg(cont, i, w_mesh, "gaussian")

    # Compute spectral dispersions using different resolution functions.
    disp_rect = spectral_disp(cont, i, w_mesh, avg_rect, "rectangle")
    disp_lorentz = spectral_disp(cont, i, w_mesh, avg_lorentz, "lorentzian")
    disp_gauss = spectral_disp(cont, i, w_mesh, avg_gauss, "gaussian")

    # Compute two-point correlators using different resolution functions.
    corr_rect = spectral_corr(cont, i, w_mesh, avg_rect, "rectangle")
    corr_lorentz = spectral_corr(cont, i, w_mesh, avg_lorentz, "lorentzian")
    corr_gauss = spectral_corr(cont, i, w_mesh, avg_gauss, "gaussian")

    # On master node, save results to an archive
    if mpi.is_master_node():
        with HDFArchive("results.h5", 'a') as ar:
            ar.create_group(str(i))
            gr = ar[str(i)]
            gr['w_mesh'] = w_mesh
            gr['avg_rect'] = avg_rect
            gr['avg_lorentz'] = avg_lorentz
            gr['avg_gauss'] = avg_gauss
            gr['disp_rect'] = disp_rect
            gr['disp_lorentz'] = disp_lorentz
            gr['disp_gauss'] = disp_gauss
            gr['corr_rect'] = corr_rect
            gr['corr_lorentz'] = corr_lorentz
            gr['corr_gauss'] = corr_gauss
