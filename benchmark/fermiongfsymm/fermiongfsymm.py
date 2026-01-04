##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2026 Igor Krivenko
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

#
# Comparison of results for the same symmetric model spectrum obtained via
# FermionGf and FermionGfSymm kernels.
#

from h5 import HDFArchive
from triqs.gf import GfImTime, GfReFreq
from triqs.gf.descriptors import Function
import triqs.utility.mpi as mpi
from som import Som, fill_refreq, reconstruct, compute_tail
from som.version import som_hash
from scipy.integrate import quad
import numpy as np


def print_master(msg):
    if mpi.rank == 0:
        print(msg)
    mpi.barrier(poll_msec=0)


#
# Parameters
#

beta = 10.0

n_tau = 201
n_iw = 201
n_w = 1001
tail_max_order = 10
energy_window = (-5, 5)


def g_tau_model_symm(tau):
    def dos(e): return np.sqrt(4 - e**2) / (2 * np.pi)
    def kern(e): return -np.exp(-tau * e) / (1 + np.exp(-beta*e))
    return quad(lambda e: dos(e) * kern(e), -2, 2, points=[0])[0]


norms = [1.0]

abs_error = 1e-4
np.random.seed(1000)

run_params = {}
run_params['verbosity'] = 2
run_params['adjust_l'] = False
run_params['t'] = 100
run_params['f'] = 100
run_params['l'] = 2000
run_params['make_histograms'] = False

if mpi.is_master_node():
    arch = HDFArchive('fermiongfsymm.h5', 'w')
    arch["som_git_hash"] = som_hash
    arch["run_params"] = run_params


def run_som_and_save(kind, g, error_bars, norms, energy_window):
    cont = Som(g, error_bars, kind=kind, norms=norms)
    cont.accumulate(energy_window=energy_window, **run_params)
    cont.compute_final_solution()

    g_w = GfReFreq(window=energy_window, n_points=n_w, indices=[0])
    fill_refreq(g_w, cont)

    g_rec = g.copy()
    reconstruct(g_rec, cont)

    tail = compute_tail(tail_max_order, cont)

    if mpi.is_master_node():
        arch.create_group(kind)
        gr = arch[kind]
        gr["input"] = g
        gr["output"] = g_w
        gr["rec"] = g_rec
        gr["tail"] = tail
        gr["solutions"] = cont.solutions


g_symm_tau = GfImTime(beta=beta, n_points=n_tau, indices=[0])

g_symm_tau << Function(g_tau_model_symm)
g_symm_tau.data[:] += abs_error * 2 * (np.random.rand(n_tau, 1, 1) - 0.5)

error_bars_tau = g_symm_tau.copy()
error_bars_tau.data[:] = abs_error

print_master("================")
print_master("FermionGf kernel")
print_master("================")

run_som_and_save("FermionGf", g_symm_tau, error_bars_tau, norms, (-5, 5))

print_master("====================")
print_master("FermionGfSymm kernel")
print_master("====================")

run_som_and_save("FermionGfSymm", g_symm_tau, error_bars_tau, norms, (-5, 5))
