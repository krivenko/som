import numpy as np
from h5 import HDFArchive
from triqs.gf import *
from triqs.gf.descriptors import *
import triqs.utility.mpi as mpi
from som import Som, fill_refreq, compute_tail, reconstruct
import time

# Parameters
beta = 20
indices = [0]
D = 1.0

n_iw = 200
n_w = 1000
tail_max_order = 9

abs_error = 0.001
energy_window = (-5, 5)
refreq_mesh = MeshReFreq(*energy_window, n_w)

accumulate_params = {'energy_window' : energy_window}
accumulate_params['verbosity'] = 2
accumulate_params['adjust_l'] = False
accumulate_params['t'] = 200
accumulate_params['f'] = 200
accumulate_params['l'] = 1000
accumulate_params['make_histograms'] = True

good_chi_rel = 2.0
good_chi_abs = np.inf

cfs_cc_params = {'refreq_mesh' : refreq_mesh,
                 'good_chi_rel' : good_chi_rel,
                 'good_chi_abs' : good_chi_abs}
cfs_cc_params['verbosity'] = 2
cfs_cc_params['max_iter'] = 10
cfs_cc_params['unity_sum_coeff'] = 1e6
cfs_cc_params['amp_penalty_max'] = 1e6
cfs_cc_params['amp_penalty_divisor'] = 100
cfs_cc_params['der_penalty_init'] = 1e-3
cfs_cc_params['der_penalty_coeff'] = 2.0

# SOCC default model
w = np.array([float(w) for w in refreq_mesh])
default_model = np.exp(-(w ** 2) / (2 * (D/2)**2)) / np.sqrt(2*np.pi*(D/2)**2)

cfs_cc_params['default_model'] = default_model
cfs_cc_params['default_model_weights'] = 1e-4 * np.ones(n_w)

print("--- Prepare input ---")
g_iw  = GfImFreq(beta = beta, n_points = n_iw, indices = indices)
g_iw << SemiCircular(D)

prng = np.random.RandomState(123456789)
g_iw.data[:] += abs_error * 2*(prng.rand(*g_iw.data.shape) - 0.5)
g_iw.data[:] = 0.5*(g_iw.data[:,:,:] + np.conj(g_iw.data[::-1,:,:]))

S_iw = g_iw.copy()
S_iw.data[:] = abs_error

def make_output(cont):
    g_iw_rec = g_iw.copy()
    reconstruct(g_iw_rec, cont)

    g_w = GfReFreq(window = energy_window, n_points = n_w, indices = indices)
    fill_refreq(g_w, cont)

    g_tail = compute_tail(tail_max_order, cont)

    return {'g_iw_rec' : g_iw_rec, 'g_w' : g_w, 'g_tail' : g_tail}

print("--- Accumulate particular solutions ---")
accumulate_time = time.perf_counter()
cont = Som(g_iw, S_iw)
cont.accumulate(**accumulate_params)
accumulate_time = time.perf_counter() - accumulate_time

print("--- Run compute_final_solution() ---")
cfs_time = time.perf_counter()
chi2 = cont.compute_final_solution(good_chi_rel=good_chi_rel,
                                   good_chi_abs=good_chi_abs)
cfs_time = time.perf_counter() - cfs_time

print("--- Generate SOM output ---")
som_output = make_output(cont)

print("--- Run compute_final_solution_cc() ---")
cfs_cc_time = time.perf_counter()
chi2_cc = cont.compute_final_solution_cc(**cfs_cc_params)
cfs_cc_time = time.perf_counter() - cfs_cc_time

print("--- Generate SOCC output ---")
socc_output = make_output(cont)

print("--- Save results ---")
if mpi.is_master_node():
    with HDFArchive('consistent_constraints.h5', 'w') as arch:
        arch.create_group('input')
        # Save input data
        input_gr = arch['input']
        input_gr['abs_error'] = abs_error
        input_gr['D'] = D
        input_gr['g_iw'] = g_iw
        input_gr['S_iw'] = S_iw
        # accumulate()
        arch['accumulate_params'] = cont.last_accumulate_parameters
        arch['accumulate_time'] = accumulate_time
        arch['histograms'] = cont.histograms
        # SOM output
        arch.create_group('som_output')
        som_output_gr = arch['som_output']
        som_output_gr['params'] = {'good_chi_rel' : good_chi_rel,
                                   'good_chi_abs' : good_chi_abs}
        som_output_gr['chi2'] = chi2
        som_output_gr['elapsed_time'] = cfs_time
        som_output_gr['g_iw_rec'] = som_output['g_iw_rec']
        som_output_gr['g_w'] = som_output['g_w']
        som_output_gr['g_tail'] = som_output['g_tail']
        # SOCC output
        arch.create_group('socc_output')
        socc_output_gr = arch['socc_output']
        socc_output_gr['params'] = cfs_cc_params
        socc_output_gr['chi2'] = chi2_cc
        socc_output_gr['elapsed_time'] = cfs_cc_time
        socc_output_gr['g_iw_rec'] = socc_output['g_iw_rec']
        socc_output_gr['g_w'] = socc_output['g_w']
        socc_output_gr['g_tail'] = socc_output['g_tail']
