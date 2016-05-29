# Reconstruction of a semicircular DOS

# Import some TRIQS modules
from pytriqs.gf.local import *
from pytriqs.gf.local.descriptors import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi

# Import main SOM class
from triqs_som.som import Som

beta = 20                   # Inverse temperature
D = 1.0                     # Half-bandwidth

n_iw = 200                  # Number of Matsubara frequencies
n_tau = 500                 # Number of tau-slices

n_w = 801                   # Number of energy slices for the solution
energy_window = (-3.0,3.0)  # Energy window to search the solution in

# Prepare input data
g_iw  = GfImFreq(beta = beta, n_points = n_iw, indices = [0])
g_iw << SemiCircular(D)

# We use the imaginary-time representation for continuation
g_tau = GfImTime(beta = beta, n_points = n_tau, indices = [0])
g_tau << InverseFourier(g_iw)

# Add some noise to the input
from numpy.random import rand, seed
seed(37283) # Seed RNG with the same value on all MPI ranks

noise = 1e-4
g_tau.data[:] += noise * (2*rand(*g_tau.data.shape) - 1)

# Set the weight function S to a constant (all points of g_tau are equally important)
S_tau = g_tau.copy()
S_tau.data[:] = 1.0

# Construct a SOM object
cont = Som(g_tau, S_tau, kind = "FermionGf")

# run() parameters
run_params = {'energy_window' : energy_window}
# Verbosity level
run_params['verbosity'] = 3
# Do not adjust the number of global updates
run_params['adjust_f'] = False
# Do not adjust the number of particular solutions to be accumulated
run_params['adjust_l'] = False
# Number of local updates per global update
run_params['t'] = 500
# Starting number of global updates (can be increased by the adjustment procedure)
run_params['f'] = 100
# Number of particular solutions
run_params['l'] = 1000
# Maximum number of rectangles to represent a solution
run_params['max_rects'] = 100
# Accumulate histogram of the objective function values
run_params['make_histograms'] = True

# Run!
cont.run(**run_params)

# Evaluate the solution on an energy mesh
# NB: we can use *any* energy window at this point, not necessarily that from run_params
g_w = GfReFreq(window = (-4.0,4.0), n_points = n_w, indices = [0])
g_w << cont

# G(\tau) reconstructed from the solution
g_rec_tau = g_tau.copy()
g_rec_tau << cont

# On master node, save results to an archive
if mpi.is_master_node():
    with HDFArchive("example.h5",'w') as ar:
        ar['g_tau'] = g_tau
        ar['g_rec_tau'] = g_rec_tau
        ar['g_w'] = g_w
