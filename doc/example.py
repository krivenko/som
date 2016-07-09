# Import some TRIQS modules
from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi

# Import main SOM class
from triqs_som.som import Som

n_tau = 500                 # Number of tau-slices for the input GF

n_w = 801                   # Number of energy slices for the solution
energy_window = (-4.0,4.0)  # Energy window to search the solution in


# Parameters for Som.run()
run_params = {'energy_window' : energy_window}
# Verbosity level
run_params['verbosity'] = 3
# Adjust the number of global updates
run_params['adjust_f'] = True
# Adjust the number of particular solutions to be accumulated
run_params['adjust_l'] = True
# Number of local updates per global update
run_params['t'] = 500
# Accumulate histogram of the objective function values
run_params['make_histograms'] = True

# Read G_tau from archive
g_tau = HDFArchive('example.h5', 'r')['g_tau']

# Prepare input data: reduce the number of \tau-slices
g_input = rebinning_tau(g_tau, n_tau)

# Set the weight function S to a constant (all points of g_tau are equally important)
S = g_input.copy()
S.data[:] = 1.0

# Construct a SOM object
cont = Som(g_input, S, kind = "FermionGf")

# Run!
cont.run(**run_params)

# Evaluate the solution on an energy mesh
# NB: we can use *any* energy window at this point, not necessarily that from run_params
g_w = GfReFreq(window = (-5.0,5.0), n_points = n_w, indices = [0,1])
g_w << cont

# G(\tau) reconstructed from the solution
g_rec_tau = g_input.copy()
g_rec_tau << cont

# On master node, save results to an archive
if mpi.is_master_node():
    with HDFArchive("results.h5",'w') as ar:
        ar['g_tau'] = g_tau
        ar['g_rec_tau'] = g_rec_tau
        ar['g_w'] = g_w
        ar['histograms'] = cont.histograms
