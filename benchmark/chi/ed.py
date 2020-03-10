from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.gf.local.descriptors import *
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import AtomDiag
from pytriqs.applications.impurity_solvers.cthyb.cthyb import *
from scipy.special import sph_in
from itertools import product
import numpy as np

from params import *

n_w = 1001
w_window = (-2.5,2.5)
w_offset = 1e-3

spin_names=('up','dn')

# Hamiltonian
H = U*n('up',0)*n('dn',0)
H += ed*(n('up',0) + n('dn',0))
H += h*(n('up',0) - n('dn',0))
for sn in spin_names:
    for i, (e, v) in enumerate(zip(eps,V)):
        bath_i = i+1
        H += e*n(sn,bath_i)
        H += v*c_dag(sn,0)*c(sn,bath_i) + v*c_dag(sn,bath_i)*c(sn,0)

# Fundamental operator set
fops = set()
for monomial in H:
    for ind in monomial[0]:
        fops.add(tuple(ind[1]))

# Diagonalize

ed = AtomDiag(H, list(fops))

print "Found", ed.n_blocks, "invariant subspaces"

Z = partition_function(ed, beta)
energies = ed.energies

# GF computation
G_tau = BlockGf(name_block_generator =
                [(sn,GfImTime(beta = beta, indices = [0], n_points = n_tau)) for sn in spin_names])
G_iw  = BlockGf(name_block_generator =
                [(sn,GfImFreq(beta = beta, indices = [0], n_points = n_iw)) for sn in spin_names])
G_l   = BlockGf(name_block_generator =
                [(sn,GfLegendre(beta = beta, indices = [0], n_points = n_l)) for sn in spin_names])
G_w   = BlockGf(name_block_generator =
                [(sn,GfReFreq(window = w_window, indices = [0], n_points = n_w)) for sn in spin_names])

for sn in spin_names:
    print "Computing GF_%s_%s" %(sn,sn)
    g_tau = G_tau[sn]
    g_iw  = G_iw[sn]
    g_l   = G_l[sn]
    g_w   = G_w[sn]
    g_terms = []
    op_linear_index = ed.fops.index([sn,0])
    for i_block in range(ed.n_blocks):
        j_block = ed.cdag_connection(op_linear_index, i_block)
        if j_block == -1: continue

        c_dag_mat = ed.cdag_matrix(op_linear_index, i_block)
        c_mat = ed.c_matrix(op_linear_index, j_block)

        for i in range(c_mat.shape[0]):
            E_i = energies[i_block][i]
            w_i = np.exp(-beta*E_i) / Z
            for j in range(c_mat.shape[1]):
                E_j = energies[j_block][j]
                w_j = np.exp(-beta*E_j) / Z

                res = (w_i + w_j)*c_mat[i,j]*c_dag_mat[j,i]
                if(abs(res) > np.finfo(float).eps):
                    dE = E_j - E_i
                    g_terms.append((dE,res))

    for (dE, res) in g_terms:
        # g_iw
        g_iw << g_iw + res*inverse(iOmega_n - dE)

        # g_tau
        fermi_fact = -1/(1 + (np.exp(-dE*beta) if dE>0 else np.exp(dE*beta)))
        for tau_n, tau in enumerate(g_tau.mesh):
            g_tau.data[tau_n,:,:] += res * fermi_fact * (np.exp(-dE*tau) if dE>0 else np.exp(dE*(beta-tau)))

        # g_l
        x = beta*dE/2
        coeff = -beta / (2 * np.cosh(x))
        sph_in_array = sph_in(len(g_l.mesh)-1,abs(x))[0]
        for l in g_l.mesh:
            l = int(l.real)
            g_l.data[l,:,:] += res * coeff*np.sqrt(2*l+1)*((-np.sign(x))**l)*sph_in_array[l]

        # g_w
        g_w << g_w + res * inverse(Omega - dE + w_offset*1j)

    # Set tails
    for g in (g_tau, g_iw, g_w):
        g.tail.zero()
        g.tail[1] = np.array([[sum(map(lambda t: t[1], g_terms))]])
        g.tail[2] = np.array([[sum(map(lambda t: t[1]*t[0], g_terms))]])
        g.tail[3] = np.array([[sum(map(lambda t: t[1]*t[0]*t[0], g_terms))]])
        g.tail.mask.fill(3)

# Chi computation
chi_tau = GfImTime(beta = beta, indices = [0,1], n_points = n_tau, statistic = "Boson")
chi_iw  = GfImFreq(beta = beta, indices = [0,1], n_points = n_iw, statistic = "Boson")
chi_l   = GfLegendre(beta = beta, indices = [0,1], n_points = n_l, statistic = "Boson")
chi_w   = GfReFreq(window = w_window, indices = [0,1], n_points = n_w)

for sn1, sn2 in product(spin_names,spin_names):
    print "Computing chi_%s_%s" %(sn1,sn2)
    g_terms = []
    op_index1 = ed.fops.index([sn1,0])
    op_index2 = ed.fops.index([sn2,0])
    chi_index1 = 0 if op_index1 == 0 else 1
    chi_index2 = 0 if op_index2 == 0 else 1
    for block in range(ed.n_blocks):
        if ed.c_connection(op_index1, block) == -1 or ed.c_connection(op_index2, block) == -1: continue

        c1_mat = np.matrix(ed.c_matrix(op_index1, block))
        c2_mat = np.matrix(ed.c_matrix(op_index2, block))
        n1_mat = c1_mat.getH() * c1_mat
        n2_mat = c2_mat.getH() * c2_mat

        for i in range(n1_mat.shape[0]):
            E_i = energies[block][i]
            w_i = np.exp(-beta*E_i) / Z
            for j in range(n1_mat.shape[1]):
                E_j = energies[block][j]
                w_j = np.exp(-beta*E_j) / Z

                dE = E_j - E_i
                if abs(dE) > np.finfo(float).eps:
                    res = (w_j - w_i)*n1_mat[i,j]*n2_mat[j,i]
                else:
                    dE = 0
                    res = beta * w_i*n1_mat[i,j]*n2_mat[j,i]
                if abs(res) > np.finfo(float).eps:
                    g_terms.append((dE,res))

    for (dE, res) in g_terms:
        # chi_iw
        if dE == 0:
            chi_iw.data[n_iw-1,chi_index1,chi_index2] += res
        else:
            chi_iw[chi_index1,chi_index2] << chi_iw[chi_index1,chi_index2] + res*inverse(iOmega_n - dE)

        # chi_tau
        if dE == 0:
            chi_tau.data[:,chi_index1,chi_index2] += res / beta
        else:
            bose_fact = -1/(1 - np.exp(-dE*beta)) if dE>0 else -1/(np.exp(dE*beta) - 1)
            for tau_n, tau in enumerate(g_tau.mesh):
                chi_tau.data[tau_n,chi_index1,chi_index2] += res * bose_fact * (np.exp(-dE*tau) if dE>0 else np.exp(dE*(beta-tau)))

        # chi_l
        if dE == 0:
            chi_l.data[0,chi_index1,chi_index2] += res
        else:
            x = beta*dE/2
            coeff = -beta / (2 * np.sinh(x))
            sph_in_array = sph_in(len(g_l.mesh)-1,abs(x))[0]
            for l in g_l.mesh:
                l = int(l.real)
                chi_l.data[l,chi_index1,chi_index2] += res * coeff*np.sqrt(2*l+1)*((-np.sign(x))**l)*sph_in_array[l]

        # chi_w
        if dE == 0:
            chi_w.data[(n_w-1)/2,chi_index1,chi_index2] += res
        else:
            chi_w[chi_index1,chi_index2] << chi_w[chi_index1,chi_index2] + res * inverse(Omega - dE + w_offset*1j)

    # Set tails
    for chi in (chi_tau, chi_iw, chi_w):
        chi.tail.zero()
        chi.tail[1][chi_index1,chi_index2] = sum(map(lambda t: t[1], g_terms))
        chi.tail[2][chi_index1,chi_index2] = sum(map(lambda t: t[1]*t[0], g_terms))
        chi.tail[3][chi_index1,chi_index2] = sum(map(lambda t: t[1]*t[0]*t[0], g_terms))
        chi.tail.mask.fill(3)

with HDFArchive(chi_ed_filename,'w') as arch:
    arch['energies'] = energies
    arch['Z'] = Z
    arch['G_tau'] = G_tau
    arch['G_iw'] = G_iw
    arch['G_l'] = G_l
    arch['G_w'] = G_w
    arch['chi_tau'] = chi_tau
    arch['chi_iw'] = chi_iw
    arch['chi_l'] = chi_l
    arch['chi_w'] = chi_w
