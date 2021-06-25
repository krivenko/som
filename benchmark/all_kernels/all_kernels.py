##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
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
# Run SOM simulations with all implemented kernels using synthetic input data
#

import time
from h5 import HDFArchive
from triqs.gf import *
from triqs.gf.descriptors import *
import triqs.utility.mpi as mpi
from som import Som, reconstruct
from som.version import som_hash
from scipy.integrate import quad
from scipy.special import spherical_in
import numpy as np

#
# Parameters
#

beta = 10
indices = [0, 1]

n_iw = 201
n_tau = 2001
n_l = 50
n_w = 1000
tail_max_order = 9

g_norms = np.array([0.7, 0.9])
chi_norms = np.array([1.2, 1.4])
chi_auto_norms = np.array([1.4, 1.6])
g_zt_norms = np.array([0.4, 0.6])

som_params = {}
som_params['verbosity'] = 2
som_params['adjust_l'] = False
som_params['t'] = 100
som_params['f'] = 100
som_params['l'] = 25
som_params['make_histograms'] = True

if mpi.is_master_node():
    arch = HDFArchive('all_kernels.np%d.h5' % mpi.world.size, 'w')
    arch["som_params"] = som_params
    arch["som_git_hash"] = som_hash

#
# Auxiliary functions
#

def print_master(msg):
    if mpi.rank == 0: print(msg)

def quad_complex(f, a, b, **kwargs):
    return quad(lambda x: f(x).real, a, b, **kwargs)[0] + 1j*quad(lambda x: f(x).imag, a, b, **kwargs)[0]

def dos(e, e0):
    return np.sqrt(4 - (e - e0)**2) / (2 * np.pi)

def run_som_and_save(kind, mesh, g, S, norms, energy_window):
    start_time = time.perf_counter()
    cont = Som(g, S, kind = kind, norms = norms)
    cont.accumulate(energy_window = energy_window, **som_params)
    cont.accumulate(energy_window = energy_window, **som_params)
    cont.compute_final_solution()
    g_rec = g.copy()
    reconstruct(g_rec, cont)
    g_w = GfReFreq(window = energy_window,
                   n_points = n_w,
                   indices = indices)
    cont.fill_observable(g_w)
    tail = cont.compute_tail(tail_max_order)
    elapsed_time = time.perf_counter() - start_time
    print_master(f"Elapsed time: %f s" % elapsed_time)
    if mpi.is_master_node():
        arch[kind].create_group(mesh)
        gr = arch[kind][mesh]
        gr["input"] = g
        gr["rec"] = g_rec
        gr["output"] = g_w
        gr["output_tail"] = tail
        gr["histograms"] = cont.histograms
        gr["solutions"] = cont.solutions
        gr["elapsed_time"] = elapsed_time

print_master("=================")
print_master("FermionGf kernels")
print_master("=================")

def g_iw_model(iw):
    kern = lambda e: 1 / (iw - e) / np.pi
    return np.diag(chi_norms) * quad_complex(lambda e: dos(e, 1) * kern(e), -1, 3, points = [0])
g_iw = GfImFreq(beta = beta, statistic = "Fermion", n_points = n_iw, indices = indices)
g_iw << Function(g_iw_model)

def g_tau_model(tau):
    kern = lambda e: -np.exp(-tau * e) / (1 + np.exp(-beta*e))
    return np.diag(chi_norms) * quad(lambda e: dos(e, 1) * kern(e), -1, 3, points = [0])[0].real
g_tau = GfImTime(beta = beta, statistic = "Fermion", n_points = n_tau, indices = indices)
g_tau << Function(g_tau_model)

def g_l_model(l):
    kern = lambda e: -beta * sqrt(2*l+1) *((-np.sign(e)) ** l) * \
                     spherical_in(l, np.abs(e)*beta/2) / (2 * np.cosh(e*beta/2))
    return np.diag(chi_norms) * quad(lambda e: dos(e, 1) * kern(e), -1, 3, points = [0])[0].real
g_l = GfLegendre(beta = beta, statistic = "Fermion", n_points = n_l, indices = indices)
g_l << Function(g_l_model)

if mpi.is_master_node():
    arch.create_group("FermionGf")
    arch["FermionGf"]["norms"] = g_norms

for mesh, g in (("imfreq", g_iw),
                ("imtime", g_tau),
                ("legendre", g_l)):
    print_master("-"*len(mesh))
    print_master(mesh)
    print_master("-"*len(mesh))
    S = g.copy()
    if mesh == "legendre": S.data[:] = 1.0
    run_som_and_save("FermionGf", mesh, g, S, g_norms, (-5,5))

print_master("=================")
print_master("BosonCorr kernels")
print_master("=================")

def chi_iw_model(iw):
    if iw == 0:
        kern = lambda e: 1 / np.pi
    else:
        kern = lambda e: -e / (iw - e) / np.pi
    return np.diag(chi_norms) * quad_complex(lambda e: dos(e, 1) * kern(e), -1, 3, points = [0])
chi_iw = GfImFreq(beta = beta, statistic = "Boson", n_points = n_iw, indices = indices)
chi_iw << Function(chi_iw_model)

def chi_tau_model(tau):
    kern = lambda e: e * np.exp(-tau * e) / (1 - np.exp(-beta*e)) / np.pi
    return np.diag(chi_norms) * quad(lambda e: dos(e, 1) * kern(e), -1, 3, points = [0])[0].real
chi_tau = GfImTime(beta = beta, statistic = "Boson", n_points = n_tau, indices = indices)
chi_tau << Function(chi_tau_model)

def chi_l_model(l):
    kern = lambda e: beta * sqrt(2*l+1) * e *((-np.sign(e)) ** l) * \
                     spherical_in(l, np.abs(e)*beta/2) / np.sinh(np.abs(e)*beta/2) / (2*np.pi)
    return np.diag(chi_norms) * quad(lambda e: dos(e, 1) * kern(e), -1, 3, points = [0])[0].real
chi_l = GfLegendre(beta = beta, statistic = "Boson", n_points = n_l, indices = indices)
chi_l << Function(chi_l_model)

if mpi.is_master_node():
    arch.create_group("BosonCorr")
    arch["BosonCorr"]["norms"] = chi_norms

for mesh, chi in (("imfreq", chi_iw),
                  ("imtime", chi_tau),
                  ("legendre", chi_l)):
    print_master("-"*len(mesh))
    print_master(mesh)
    print_master("-"*len(mesh))
    S = chi.copy()
    if mesh == "legendre": S.data[:] = 1.0
    run_som_and_save("BosonCorr", mesh, chi, S, chi_norms, (-5,5))

print_master("=====================")
print_master("BosonAutoCorr kernels")
print_master("=====================")

def chi_auto_iw_model(iw):
    if iw == 0:
        kern = lambda e: 2 / np.pi
    else:
        kern = lambda e: 2 * (e**2) / (e**2 + iw.imag**2) / np.pi
    return np.diag(chi_auto_norms) * quad_complex(lambda e: 2*dos(e, 0) * kern(e), 0, 2)
chi_auto_iw = GfImFreq(beta = beta, statistic = "Boson", n_points = n_iw, indices = indices)
chi_auto_iw << Function(chi_auto_iw_model)

def chi_auto_tau_model(tau):
    kern = lambda e: e * (np.exp(-tau * e) + np.exp(-(beta-tau) * e)) / (1 - np.exp(-beta*e)) / np.pi
    return np.diag(chi_auto_norms) * quad(lambda e: 2*dos(e, 0) * kern(e), 0, 2)[0].real
chi_auto_tau = GfImTime(beta = beta, statistic = "Boson", n_points = n_tau, indices = indices)
chi_auto_tau << Function(chi_auto_tau_model)

def chi_auto_l_model(l):
    kern = lambda e: beta * (1 + (-1) ** l) * sqrt(2*l+1) * e * spherical_in(l, e*beta/2) / np.sinh(e*beta/2) / (2*np.pi)
    return np.diag(chi_auto_norms) * quad(lambda e: 2*dos(e, 0) * kern(e), 0, 2)[0].real
chi_auto_l = GfLegendre(beta = beta, statistic = "Boson", n_points = n_l, indices = indices)
chi_auto_l << Function(chi_auto_l_model)

if mpi.is_master_node():
    arch.create_group("BosonAutoCorr")
    arch["BosonAutoCorr"]["norms"] = chi_auto_norms

for mesh, chi_auto in (("imfreq", chi_auto_iw),
                       ("imtime", chi_auto_tau),
                       ("legendre", chi_auto_l)):
    print_master("-"*len(mesh))
    print_master(mesh)
    print_master("-"*len(mesh))
    S = chi_auto.copy()
    if mesh == "legendre": S.data[:] = 1.0
    run_som_and_save("BosonAutoCorr",
                     mesh,
                     chi_auto,
                     S,
                     chi_auto_norms,
                     (0,5))

print_master("================")
print_master("ZeroTemp kernels")
print_master("================")

def g_zt_iw_model(iw):
    kern = lambda e: 1 / (iw - e)
    return np.diag(g_zt_norms) * quad_complex(lambda e: 2*dos(e, 0) * kern(e), 0, 2)
g_zt_iw = GfImFreq(beta = beta, n_points = n_iw, indices = indices)
g_zt_iw << Function(g_zt_iw_model)

def g_zt_tau_model(tau):
    kern = lambda e: -np.exp(-tau * e)
    return np.diag(g_zt_norms) * quad(lambda e: 2*dos(e, 0) * kern(e), 0, 2)[0].real
g_zt_tau = GfImTime(beta = beta, n_points = n_tau, indices = indices)
g_zt_tau << Function(g_zt_tau_model)

def g_zt_l_model(l):
    kern = lambda e: beta * ((-1) ** (l+1)) * sqrt(2*l+1) * spherical_in(l, e*beta/2)*np.exp(-e*beta/2)
    return np.diag(g_zt_norms) * quad(lambda e: 2*dos(e, 0) * kern(e), 0, 2)[0].real
g_zt_l = GfLegendre(beta = beta, n_points = n_l, indices = indices)
g_zt_l << Function(g_zt_l_model)

if mpi.is_master_node():
    arch.create_group("ZeroTemp")
    arch["ZeroTemp"]["norms"] = g_zt_norms

for mesh, g_zt in (("imfreq", g_zt_iw),
                   ("imtime", g_zt_tau),
                   ("legendre", g_zt_l)):
    print_master("-"*len(mesh))
    print_master(mesh)
    print_master("-"*len(mesh))
    S = g_zt.copy()
    if mesh == "legendre": S.data[:] = 1.0
    run_som_and_save("ZeroTemp", mesh, g_zt, S, g_zt_norms, (0,5))
