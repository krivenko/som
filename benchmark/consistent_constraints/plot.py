from h5 import HDFArchive
from triqs.gf import *
from triqs.gf.descriptors import *
from triqs.plot.mpl_interface import oplot
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from sys import argv

arch_name = argv[1]

arch = HDFArchive(arch_name, 'r')
pp = PdfPages('consistent_constraints.pdf')

D = arch['input']['D']
energy_window = arch['accumulate_params']['energy_window']

g_w_som = arch['som_output']['g_w']
chi_som = arch['som_output']['chi2'][0]

g_w_socc = arch['socc_output']['g_w']
chi_socc = arch['socc_output']['chi2'][0]

socc_iterations = arch['socc_output']['iterations']
N_iter = len(socc_iterations)
N_k = len(socc_iterations[0]['Q_k'])

A_k, Q_k = np.zeros((N_k, N_iter)), np.zeros((N_k, N_iter))
Ap_k, D_k = np.zeros((N_k - 1, N_iter)), np.zeros((N_k - 1, N_iter))
App_k, B_k = np.zeros((N_k - 2, N_iter)), np.zeros((N_k - 2, N_iter))
for i, it in enumerate(socc_iterations):
    A_k[:, i] = it['A_k']
    Q_k[:, i] = it['Q_k']
    Ap_k[:, i] = it['Ap_k']
    D_k[:, i] = it['D_k']
    App_k[:, i] = it['App_k']
    B_k[:, i] = it['B_k']

g_w_ref = g_w_som.copy()
g_w_ref << SemiCircular(D)

w = np.array([float(w) for w in g_w_ref.mesh])
default_model = arch['socc_output']['params']['default_model']

fig = plt.figure()
fig.suptitle("Spectral functions")

# Plot referenced spectral function
oplot(g_w_ref, mode = 'S', label = "Reference")

# Plot default model
plt.plot(w, default_model, label = "SOCC default model")

# Plot SOM spectral function
oplot(g_w_som, mode = 'S', label = f"SOM, $\\chi^2=%f$" % chi_som)

# Plot SOCC spectral function
oplot(g_w_socc, mode = 'S', label = f"SOCC, $\\chi^2=%f$" % chi_socc)

plt.xlim(energy_window)
plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")

pp.savefig(fig)

#
# Plot SOCC regularization parameters
#

# Amplitude regularization parameters

fig, axes = plt.subplots(1, 2, sharey = 'row')
fig.suptitle("Convergence of SOCC amplitude regularization parameters")

A_k_image = axes[0].imshow(A_k, extent = (0, N_iter, w[-1], w[0]))
axes[0].set_title("$A_k$")
axes[0].set_xlabel("# iteration")
axes[0].set_ylabel(r"$\omega_k$")
plt.colorbar(A_k_image, ax=axes[0], orientation='horizontal')

Q_k_image = axes[1].imshow(Q_k, extent = (0, N_iter, w[-1], w[0]))
axes[1].set_title("$Q_k$")
axes[1].set_xlabel("# iteration")
axes[1].set_ylabel(r"$\omega_k$")
plt.colorbar(Q_k_image, ax=axes[1],
             orientation='horizontal',
             norm=LogNorm(vmin=Q_k.min(), vmax=Q_k.max()))

pp.savefig(fig)

# Derivative regularization parameters

fig, axes = plt.subplots(1, 2, sharey = 'row')
fig.suptitle("Convergence of SOCC derivative regularization parameters")

Ap_k_image = axes[0].imshow(Ap_k, extent = (0, N_iter, w[-1], w[1]))
axes[0].set_title("$A'_k$")
axes[0].set_xlabel("# iteration")
axes[0].set_ylabel(r"$\omega_k$")
plt.colorbar(Ap_k_image, ax=axes[0], orientation='horizontal')

D_k_image = axes[1].imshow(D_k, extent = (0, N_iter, w[-1], w[1]))
axes[1].set_title("$D_k$")
axes[1].set_xlabel("# iteration")
axes[1].set_ylabel(r"$\omega_k$")
plt.colorbar(D_k_image, ax=axes[1],
             orientation='horizontal',
             norm=LogNorm(vmin=D_k.min(), vmax=D_k.max()))

pp.savefig(fig)

# Second derivative regularization parameters

fig, axes = plt.subplots(1, 2, sharey = 'row')
fig.suptitle("Convergence of SOCC second derivative regularization parameters")

App_k_image = axes[0].imshow(App_k, extent = (0, N_iter, w[-2], w[1]))
axes[0].set_title("$A''_k$")
axes[0].set_xlabel("# iteration")
axes[0].set_ylabel(r"$\omega_k$")
plt.colorbar(App_k_image, ax=axes[0], orientation='horizontal')

B_k_image = axes[1].imshow(B_k, extent = (0, N_iter, w[-2], w[1]))
axes[1].set_title("$D_k$")
axes[1].set_xlabel("# iteration")
axes[1].set_ylabel(r"$\omega_k$")
plt.colorbar(B_k_image, ax=axes[1],
             orientation='horizontal',
             norm=LogNorm(vmin=B_k.min(), vmax=B_k.max()))

pp.savefig(fig)

pp.close()
