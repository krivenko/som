##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2024 Igor Krivenko
#
# SOM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# SOM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# SOM. If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

from h5 import HDFArchive
from triqs.gf import GfImFreq, GfImTime, GfLegendre, GfReFreq
from triqs.gf.descriptors import SemiCircular, Fourier, MatsubaraToLegendre
import triqs.utility.mpi as mpi
from som import Som, fill_refreq, compute_tail, reconstruct
import numpy as np
import time

beta = 20
indices = [0]
D = 1.0

n_iw = 200
n_tau = 500
n_l = 50
n_w = 1000
tail_max_order = 9

abs_error = [0.01, 0.005, 0.002, 0.001]

run_params = {'energy_window': (-5, 5)}
run_params['verbosity'] = 2
run_params['adjust_l'] = False
run_params['t'] = 1000
run_params['f'] = 1000
run_params['l'] = 500
run_params['make_histograms'] = True

g_iw = GfImFreq(beta=beta, n_points=n_iw, indices=indices)
g_tau = GfImTime(beta=beta, n_points=n_tau, indices=indices)
g_l = GfLegendre(beta=beta, n_points=n_l, indices=indices)

g_w = GfReFreq(window=run_params['energy_window'],
               n_points=n_w,
               indices=indices)

error_bars_iw = g_iw.copy()
error_bars_tau = g_tau.copy()
error_bars_l = g_l.copy()

g_iw_rec = g_iw.copy()
g_tau_rec = g_tau.copy()
g_l_rec = g_l.copy()

if mpi.is_master_node():
    arch = HDFArchive('semicircle.h5', 'w')
    arch['abs_errors'] = abs_error
    arch['D'] = D

for s in abs_error:
    if mpi.is_master_node():
        g_iw << SemiCircular(D)

        g_tau << Fourier(g_iw)
        g_tau.data[:] = g_tau.data[:].real

        g_l << MatsubaraToLegendre(g_iw)
        g_l.data[:] = g_l.data[:].real

        g_iw.data[:] += s * 2*(np.random.rand(*g_iw.data.shape) - 0.5)
        g_iw.data[:] = 0.5*(g_iw.data[:, :, :] + np.conj(g_iw.data[::-1, :, :]))
        g_tau.data[:] += s * 2*(np.random.rand(*g_tau.data.shape) - 0.5)
        g_l.data[:] += s * 2*(np.random.rand(*g_l.data.shape) - 0.5)

    g_iw = mpi.bcast(g_iw)
    g_tau = mpi.bcast(g_tau)
    g_l = mpi.bcast(g_l)

    error_bars_iw.data[:] = 1.0
    error_bars_tau.data[:] = 1.0
    error_bars_l.data[:] = 1.0

    if mpi.is_master_node():
        gr_name = 'abs_error_%.4f' % s
        arch.create_group(gr_name)
        abs_err_gr = arch[gr_name]

    for name, g, error_bars, g_rec in (
            ('g_iw', g_iw, error_bars_iw, g_iw_rec),
            ('g_tau', g_tau, error_bars_tau, g_tau_rec),
            ('g_l',  g_l,  error_bars_l,  g_l_rec)
            ):

        start = time.perf_counter()
        cont = Som(g, error_bars)
        cont.run(**run_params)
        exec_time = time.perf_counter() - start

        reconstruct(g_rec, cont)
        fill_refreq(g_w, cont)
        g_tail = compute_tail(tail_max_order, cont)

        if mpi.is_master_node():
            abs_err_gr.create_group(name)
            gr = abs_err_gr[name]
            gr['params'] = cont.last_accumulate_parameters
            gr['exec_time'] = exec_time
            gr['g'] = g
            gr['g_w'] = g_w
            gr['g_rec'] = g_rec
            gr['g_tail'] = g_tail
            gr['histograms'] = cont.histograms
