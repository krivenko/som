from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from triqs_som.som import Som

beta = 10
indices = [1,2]

g_tau = GfImTime(beta = beta, n_points = 200, indices = indices)
s_tau = GfImTime(beta = beta, n_points = 200, indices = indices)

g_w = GfReFreq(window = (-5,5), n_points = 1000, indices = indices)

cont = Som(g_tau, s_tau)

run_params = {'energy_window' : (-5,5)}

#cont.run(**run_params)
#g_w << cont

if mpi.is_master_node():
    arch = HDFArchive('python_gf.out.h5','w')
    arch['g_w'] = g_w
