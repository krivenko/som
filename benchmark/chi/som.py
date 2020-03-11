from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from pytriqs.applications.analytical_continuation.som import Som
import numpy as np

from sys import argv
from params import *


# Fixed parameters for Som.run()
som_params['verbosity'] = 2
som_params['make_histograms'] = True

kind, mesh = argv[1:]
mesh_suffix = {'imtime':'tau','imfreq':'iw','legendre':'l'}[mesh]

arch = HDFArchive(chi_filename,'r')
arch_ed = HDFArchive(chi_ed_filename,'r')

# Read and preprocess input data
if kind == 'FermionGf':
    G = arch['G_' + mesh_suffix]
    G_w = arch_ed['G_w']
    n_w = len(G_w.mesh)
    energy_window = (G_w.mesh.omega_min, G_w.mesh.omega_max)
elif kind == 'BosonCorr' or kind == 'BosonAutoCorr':
    if mesh != 'legendre':
        chi = arch['chi_' + mesh_suffix]
    else:
        chi = GfLegendre(beta = beta,
                         statistic = "Boson",
                         indices = [0,1],
                         n_points = n_l)
        chi << MatsubaraToLegendre(arch['chi_iw'])
    chi.data[:] = chi.data.real
    chi_w = arch_ed['chi_w'].copy()
    n_w = len(chi_w.mesh)
    energy_window = (chi_w.mesh.omega_min, chi_w.mesh.omega_max)
else:
    raise RuntimeError('Unknown observable kind '+kind)

som_params['energy_window'] = energy_window

arch_som = HDFArchive(som_filename,'a')

if kind == 'FermionGf':
    G_rec = G.copy()
    G_w = arch_ed['G_w'].copy()
    histograms = {}
    for bn, g_block in G:
        # Construct a SOM object
        cont = Som(g_block, kind = kind)
        # Run!
        cont.run(**som_params)

        G_w[bn] << cont
        G_rec[bn] << cont
        histograms[bn] = cont.histograms

    if mpi.is_master_node():
        with HDFArchive(som_filename,'a') as arch_som:
            gr = 'G_'+mesh_suffix
            if gr not in arch_som: arch_som.create_group(gr)
            arch_som[gr]['G'] = G
            arch_som[gr]['G_rec'] = G_rec
            arch_som[gr]['G_w'] = G_w
            arch_som[gr]['histograms'] = histograms
else:
    chi_rec = chi.copy()
    chi_w = arch_ed['chi_w'].copy()

    chi_iw = arch['chi_iw']
    norms = np.pi * np.array([chi_iw.data[n_iw-1,0,0].real,
                              chi_iw.data[n_iw-1,1,1].real])
    if kind == "BosonAutoCorr": norms = norms / 2

    # Construct a SOM object
    cont = Som(chi, kind = kind, norms = norms)
    # Run!
    cont.run(**som_params)

    chi_w << cont
    chi_rec << cont

    if mpi.is_master_node():
        with HDFArchive(som_filename,'a') as arch_som:
            gr = {"BosonCorr" : 'chi_', "BosonAutoCorr" : 'chi_auto_'}[kind] + mesh_suffix
            if gr not in arch_som: arch_som.create_group(gr)
            arch_som[gr]['chi'] = chi
            arch_som[gr]['chi_rec'] = chi_rec
            arch_som[gr]['chi_w'] = chi_w
            arch_som[gr]['histograms'] = cont.histograms
