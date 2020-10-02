# Import some TRIQS modules and NumPy
from h5 import HDFArchive
from triqs.gf import *
import triqs.utility.mpi as mpi
import numpy

# Import main SOM class
from som import Som

n_w = 1001                  # Number of energy slices for the solution
energy_window = (-5.0,5.0)  # Energy window to search the solution in
tail_max_order = 5          # Maximum high energy expansion order to be computed

# Parameters for Som.run()
run_params = {'energy_window' : energy_window}
# Verbosity level
run_params['verbosity'] = 3
# Number of particular solutions to accumulate
run_params['l'] = 10000
# Number of global updates
run_params['f'] = 100
# Number of local updates per global update
run_params['t'] = 50
# Accumulate histogram of the objective function values
run_params['make_histograms'] = True

# Read \chi(i\omega_n) from archive
# Could be \chi(\tau) or \chi_l as well.
chi_iw = HDFArchive('example.h5', 'r')['chi_iw']

# Set the weight function S to a constant (all points of chi_iw are equally important)
S = chi_iw.copy()
S.data[:] = 1.0

# Construct a SOM object
#
# Norms of spectral functions are known to be 3.0 for both diagonal components of \chi.
# It is possible to estimate the norms as \pi\chi(i\omega_0) (estimate affected by noise)
cont = Som(chi_iw, S, kind = "BosonCorr", norms = numpy.array([3.0, 3.0]))

# Run!
# Takes 10-15 minutes on 16 cores ...
cont.run(**run_params)

# Evaluate the solution on an energy mesh
# NB: we can use *any* energy window at this point, not necessarily that from run_params
chi_w = GfReFreq(window = (-5.0,5.0), n_points = n_w, indices = [0,1])
cont.fill_observable(chi_w)

# \chi(i\omega_n) reconstructed from the solution
chi_rec_iw = chi_iw.copy()
cont.fill_observable(chi_rec_iw)

# Compute tail coefficients
tail = cont.compute_tail(tail_max_order)

# On master node, save results to an archive
if mpi.is_master_node():
    with HDFArchive("results.h5",'w') as ar:
        ar['chi_iw'] = chi_iw
        ar['chi_rec_iw'] = chi_rec_iw
        ar['chi_w'] = chi_w
        ar['tail'] = tail
        ar['histograms'] = cont.histograms
