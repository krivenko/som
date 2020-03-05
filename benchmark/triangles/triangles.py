from pytriqs.gf import *
from pytriqs.gf.descriptors import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from som import Som
import numpy as np
import time

from dos import make_g_tau

beta = 50
indices = [0]

n_tau = 600
n_w = 1000
tail_max_order = 9

abs_error = [1e-4]

run_params = {'energy_window' : (-4.0,7.0)}
run_params['verbosity'] = 2
run_params['adjust_f'] = True
run_params['adjust_l'] = True
run_params['t'] = 500
run_params['f'] = 500
run_params['l'] = 150
run_params['make_histograms'] = True

g_tau = GfImTime(beta = beta, n_points = n_tau, indices = indices)
g_w = GfReFreq(window = run_params['energy_window'], n_points = n_w, indices = indices)
S_tau = g_tau.copy()
g_tau_rec = g_tau.copy()

if mpi.is_master_node():
    arch = HDFArchive('triangles.h5','w')
    arch['abs_errors'] = abs_error

for s in abs_error:
    if mpi.is_master_node():
        make_g_tau(g_tau)
        g_tau.data[:] += s * 2*(np.random.rand(*g_tau.data.shape) - 0.5)

    g_tau = mpi.bcast(g_tau)
    S_tau.data[:] = 1.0

    if mpi.is_master_node():
        gr_name = 'abs_error_%.4f' % s
        arch.create_group(gr_name)
        abs_err_gr = arch[gr_name]

    for name, g, S, g_rec in (('g_tau',g_tau,S_tau,g_tau_rec),):

        start = time.clock()
        cont = Som(g, S)
        cont.run(**run_params)
        exec_time = time.clock() - start

        cont.fill_observable(g_rec)
        cont.fill_observable(g_w)
        g_tail = cont.compute_tail(tail_max_order)

        if mpi.is_master_node():
            abs_err_gr.create_group(name)
            gr = abs_err_gr[name]
            gr['params'] = cont.last_run_parameters
            gr['exec_time'] = exec_time
            gr['g'] = g
            gr['g_w'] = g_w
            gr['g_rec'] = g_rec
            gr['g_tail'] = g_tail
            gr['histograms'] = cont.histograms
