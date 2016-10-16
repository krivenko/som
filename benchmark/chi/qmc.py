import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.operators import *
from cthyb_segment import Solver

from params import *

p = {}
p['n_cycles'] = 4000000
p['n_warmup_cycles'] = 1000000
p['length_cycle'] = 50
p['verbosity'] = 3
p['move_move'] = True
p['measure_gt'] = True
p['measure_gw'] = True
p['measure_gl'] = True
p['measure_nnt'] = True
p['measure_nnw'] = True

# Interaction Hamiltonian
h_int = U*n('up',0)*n('dn',0)

S = Solver(beta = beta,
           gf_struct = gf_struct,
           n_tau_g = n_tau,
           n_w = n_iw,
           n_legendre_g = n_l,
           n_tau_nn = n_tau,
           n_w_b_nn = n_iw)

delta_iw = GfImFreq(beta = beta, indices = [0], n_points = n_iw)
for e, v in zip(eps,V):
    delta_iw << delta_iw + v**2 * inverse(iOmega_n - e)
for bn, g0 in S.G0_iw:
    g0 << inverse(iOmega_n - ed - {'up':h,'dn':-h}[bn] - delta_iw)

S.solve(h_int = h_int, **p)

if mpi.is_master_node():
    with HDFArchive(chi_filename,'w') as arch:
        arch['params'] = p
        arch['delta_iw']  = delta_iw
        arch['Delta_tau'] = S.Delta_tau
        arch['G0_iw']     = S.G0_iw
        arch['G_tau']     = S.G_tau
        arch['G_iw']      = S.G_iw
        arch['G_l']       = S.G_l
        arch['chi_iw']    = S.nn_iw
        arch['chi_tau']   = S.nn_tau
