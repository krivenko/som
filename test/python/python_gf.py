from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from pytriqs.applications.analytic_continuation.som import SomCore

beta = 10
indices = [1,2]

g_tau = GfImTime(beta = beta, n_points = 200, indices = indices)
s_tau = GfImTime(beta = beta, n_points = 200, indices = indices)

g_w = GfReFreq(window = (-5,5), n_points = 1000, indices = indices)

cont = SomCore(g_tau, s_tau)

run_params = {}

cont.run(**run_params)
cont(g_w)

if mpi.is_master_node():
    arch = HDFArchive('python_gf.out.h5','w')
    arch['g_w'] = g_w
