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
# BosonCorr and BosonAutoCorr kernels.
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


def chi_tau_model_symm(tau):
    def dos(e):
        return np.sqrt(4 - e**2) / (2 * np.pi)

    def kern(e):
        if e == 0:
            return 1 / (beta * np.pi)
        else:
            return e * np.exp(-tau * e) / (1 - np.exp(-beta*e)) / np.pi
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
    arch = HDFArchive('bosonautocorr.h5', 'w')
    arch["som_git_hash"] = som_hash
    arch["run_params"] = run_params


def run_som_and_save(kind, chi, error_bars, norms, energy_window):
    cont = Som(chi, error_bars, kind=kind, norms=norms)
    cont.accumulate(energy_window=energy_window, **run_params)
    cont.compute_final_solution()

    chi_w = GfReFreq(window=energy_window, n_points=n_w, indices=[0])
    fill_refreq(chi_w, cont)

    chi_rec = chi.copy()
    reconstruct(chi_rec, cont)

    tail = compute_tail(tail_max_order, cont)

    if mpi.is_master_node():
        arch.create_group(kind)
        gr = arch[kind]
        gr["input"] = chi
        gr["output"] = chi_w
        gr["rec"] = chi_rec
        gr["tail"] = tail
        gr["solutions"] = cont.solutions


chi_symm_tau = GfImTime(beta=beta,
                        statistic="Boson",
                        n_points=n_tau,
                        indices=[0])

chi_symm_tau << Function(chi_tau_model_symm)
chi_symm_tau.data[:] += abs_error * 2 * (np.random.rand(n_tau, 1, 1) - 0.5)

error_bars_tau = chi_symm_tau.copy()
error_bars_tau.data[:] = abs_error

print_master("================")
print_master("BosonCorr kernel")
print_master("================")

run_som_and_save("BosonCorr",
                 chi_symm_tau,
                 error_bars_tau,
                 norms,
                 energy_window)

print_master("====================")
print_master("BosonAutoCorr kernel")
print_master("====================")

run_som_and_save("BosonAutoCorr",
                 chi_symm_tau,
                 error_bars_tau,
                 norms,
                 energy_window)
