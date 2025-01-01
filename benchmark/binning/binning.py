##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2025 Igor Krivenko
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
# Projection of observables on a real frequency grid with and without binning
#

import time
from h5 import HDFArchive
from triqs.gf import GfImFreq, GfReFreq
from triqs.gf.descriptors import Function
import triqs.utility.mpi as mpi
from som import Som, fill_refreq
from som.version import som_hash
from scipy.integrate import quad
import numpy as np


def print_master(msg):
    if mpi.rank == 0:
        print(msg)
    mpi.barrier(poll_msec=0)


def quad_complex(f, a, b, **kwargs):
    return quad(lambda x: f(x).real, a, b, **kwargs)[0] + \
        1j*quad(lambda x: f(x).imag, a, b, **kwargs)[0]


#
# Parameters
#

beta = 10.0

n_iw = 201
n_w_nobinning = 10001
n_w = 101


def dos(e, e0):
    return np.sqrt(4 - (e - e0)**2) / (2 * np.pi)


run_params = {}
run_params['verbosity'] = 2
run_params['adjust_l'] = False
run_params['t'] = 100
run_params['f'] = 100
run_params['l'] = 1000
run_params['make_histograms'] = False

if mpi.is_master_node():
    arch = HDFArchive('binning.h5', 'w')
    arch["som_git_hash"] = som_hash
    arch["run_params"] = run_params


def run_som_and_save(kind, g, error_bars, norms, energy_window):
    cont = Som(g, error_bars, kind=kind, norms=norms)
    cont.accumulate(energy_window=energy_window, **run_params)
    cont.compute_final_solution()

    g_w = GfReFreq(window=energy_window, n_points=n_w, indices=[0])
    g_w_nobinning = GfReFreq(window=energy_window,
                             n_points=n_w_nobinning,
                             indices=[0])

    start_time = time.perf_counter()
    fill_refreq(g_w_nobinning, cont, with_binning=False)
    nobinning_time = time.perf_counter()
    fill_refreq(g_w, cont, with_binning=True)
    binning_time = time.perf_counter()

    elapsed_time_nobinning = nobinning_time - start_time
    elapsed_time_binning = binning_time - nobinning_time

    print_master("Elapsed time: %f s without binning, %f s with binning" %
                 (elapsed_time_nobinning, elapsed_time_binning))

    if mpi.is_master_node():
        arch.create_group(kind)
        gr = arch[kind]
        gr["input"] = g_iw
        gr["output"] = g_w
        gr["output_nobinning"] = g_w_nobinning
        gr["solutions"] = cont.solutions
        gr["elapsed_time_nobinning"] = elapsed_time_nobinning
        gr["elapsed_time_binning"] = elapsed_time_binning


print_master("================")
print_master("FermionGf kernel")
print_master("================")


def g_iw_model(iw):
    def kern(e):
        return 1 / (iw - e)
    return quad_complex(lambda e: dos(e, 1) * kern(e), -1, 3, points=[0])


g_iw = GfImFreq(beta=beta, n_points=n_iw, indices=[0])
g_iw << Function(g_iw_model)
error_bars_iw = g_iw.copy()
error_bars_iw.data[:] = np.abs(error_bars_iw.data[:])

run_som_and_save("FermionGf", g_iw, error_bars_iw, [1.0], (-5, 5))

print_master("====================")
print_master("FermionGfSymm kernel")
print_master("====================")


def g_iw_model(iw):
    def kern(e):
        return 1 / (iw - e)
    return quad_complex(lambda e: dos(e, 0) * kern(e), -2, 2, points=[0])


g_iw = GfImFreq(beta=beta, n_points=n_iw, indices=[0])
g_iw << Function(g_iw_model)
error_bars_iw = g_iw.copy()
error_bars_iw.data[:] = np.abs(error_bars_iw.data[:])

run_som_and_save("FermionGfSymm", g_iw, error_bars_iw, [1.0], (-5, 5))

print_master("================")
print_master("BosonCorr kernel")
print_master("================")

chi_norms = np.array([1.2])


def chi_iw_model(iw):
    if iw == 0:
        def kern(e):
            return 1 / np.pi
    else:
        def kern(e):
            return -e / (iw - e) / np.pi
    return chi_norms * quad_complex(lambda e: dos(e, 1) * kern(e),
                                    -1, 3,
                                    points=[0])


chi_iw = GfImFreq(beta=beta,
                  statistic="Boson",
                  n_points=n_iw,
                  indices=[0])
chi_iw << Function(chi_iw_model)
error_bars_iw = chi_iw.copy()
error_bars_iw.data[:] = np.abs(error_bars_iw.data[:])

run_som_and_save("BosonCorr", chi_iw, error_bars_iw, chi_norms, (-5, 5))

print_master("====================")
print_master("BosonAutoCorr kernel")
print_master("====================")

chi_norms = np.array([1.2])


def chi_iw_model(iw):
    if iw == 0:
        def kern(e):
            return 2 / np.pi
    else:
        def kern(e):
            return 2 * (e**2) / (e**2 + iw.imag**2) / np.pi
    return chi_norms * quad_complex(lambda e: 2*dos(e, 0) * kern(e), 0, 2)


chi_iw = GfImFreq(beta=beta,
                  statistic="Boson",
                  n_points=n_iw,
                  indices=[0])
chi_iw << Function(chi_iw_model)
error_bars_iw = chi_iw.copy()
error_bars_iw.data[:] = np.abs(error_bars_iw.data[:])

run_som_and_save("BosonAutoCorr", chi_iw, error_bars_iw, chi_norms, (-5, 5))

print_master("===============")
print_master("ZeroTemp kernel")
print_master("===============")

g_zt_norms = np.array([0.4])


def g_zt_iw_model(iw):
    def kern(e):
        return 1 / (iw - e)
    return g_zt_norms * quad_complex(lambda e: 2*dos(e, 0) * kern(e), 0, 2)


g_zt_iw = GfImFreq(beta=beta, n_points=n_iw, indices=[0])
g_zt_iw << Function(g_zt_iw_model)
error_bars_zt_iw = g_zt_iw.copy()
error_bars_zt_iw.data[:] = np.abs(error_bars_zt_iw.data[:])

run_som_and_save("ZeroTemp", g_zt_iw, error_bars_zt_iw, g_zt_norms, (0, 5))
