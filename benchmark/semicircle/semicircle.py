from pytriqs.gf.local import *
from pytriqs.gf.local.descriptors import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from triqs_som.som import Som
import numpy as np

beta = 20
indices = [1,2]
D = 1.0

n_iw = 200
n_tau = 500
abs_error = 0.005

seed = 1
numpy.random.seed(seed)

g_iw = GfImFreq(beta = beta, n_points = n_iw, indices= indices)
g_iw << SemiCircular(D)

g_tau = GfImTime(beta = beta, n_points = n_tau, indices = indices)
g_tau << InverseFourier(g_iw)

g_tau.data[:] += abs_error * (np.random.rand(*g_tau.data.shape) - 0.5)

s_tau = GfImTime(beta = beta, n_points = n_tau, indices = indices)
s_tau.data[:] = abs_error

g_w = GfReFreq(window = (-5,5), n_points = 1000, indices = indices)

cont = Som(g_tau, s_tau)

run_params = {'energy_window' : (-5,5)}
run_params['verbosity'] = 3
run_params['random_seed'] = seed
run_params['n_elementary_updates'] = 500
run_params['adjust_ngu'] = True
run_params['n_global_updates'] = 50
run_params['adjust_ngu_n_solutions'] = 10;
run_params['adjust_nsol'] = True
run_params['n_solutions'] = 10
run_params['adjust_nsol_ratio'] = 0.8
run_params['make_histograms'] = True

cont.run(**run_params)
g_w << cont

# Reconstructed G(\tau)
g_tau_rec = g_tau.copy()
g_tau_rec << cont

if mpi.is_master_node():
    arch = HDFArchive('python_gf.out.h5','w')
    arch['g_w'] = g_w
    arch['g_tau_rec'] = g_tau_rec
    arch['histograms'] = cont.histograms
