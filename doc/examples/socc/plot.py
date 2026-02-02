from h5 import HDFArchive
from triqs.gf import Gf                                             # noqa: F401
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from triqs.plot.mpl_interface import oplot
import numpy as np

# Read data from archive
ar = HDFArchive("results.h5", 'r')

# Energy mesh
energy_mesh = ar["som"]["g_w"].mesh
w_points = np.array([float(w) for w in energy_mesh])

# Plot spectral functions obtained using SOM and SOCC procedures
oplot(ar["som"]["g_w"], mode='S', lw=0.8,
      label=r'SOM, $\chi^2=%f$' % ar["som"]["chi2"][0])
oplot(ar["socc"]["g_w"], mode='S', lw=0.8,
      label=r'SOCC, $\chi^2=%f$' % ar["socc"]["chi2"][0])

# Plot the default model
plt.plot(w_points, ar["socc"]["default_model"], lw=0.8,
         label="default model")

plt.xlim((energy_mesh.w_min, energy_mesh.w_max))
plt.ylim((0, 0.4))
plt.ylabel(r"$A(\omega)$")
plt.legend()
plt.show()

# Plot input and reconstructed G(\tau)
oplot(ar["som"]["g_tau"][0, 0], mode='R', lw=0.8, label=r"$G(\tau)$")
oplot(ar["som"]["g_rec_tau"][0, 0], mode='R', lw=0.8,
      label=r"$G^\mathrm{rec}(\tau)$, SOM")
oplot(ar["socc"]["g_rec_tau"][0, 0], mode='R', lw=0.8,
      label=r"$G^\mathrm{rec}(\tau)$, SOCC")

plt.xlim((0, ar["som"]["g_tau"].mesh.beta))
plt.ylabel(r"$G(\tau)$")
plt.legend(loc="lower center")
plt.show()

#
# Plot evolution of SOCC regularization parameters with iterations
#

socc_iterations = ar['socc']['iterations']
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

# Amplitude regularization parameters
fig, axes = plt.subplots(1, 2, sharey='row', figsize=(8, 3.6))
fig.suptitle("Evolution of SOCC amplitude regularization parameters")

A_k_image = axes[0].imshow(A_k, extent=(0, N_iter, w_points[-1], w_points[0]))
axes[0].set_title(r"$A_k$")
axes[0].set_xlabel(r"# iteration")
axes[0].set_ylabel(r"$\omega_k$")
plt.colorbar(A_k_image, ax=axes[0], orientation='horizontal')

Q_k_image = axes[1].imshow(Q_k, extent=(0, N_iter, w_points[-1], w_points[0]))
axes[1].set_title(r"$Q_k$")
axes[1].set_xlabel(r"# iteration")
axes[1].set_ylabel(r"$\omega_k$")
plt.colorbar(Q_k_image,
             ax=axes[1],
             orientation='horizontal',
             norm=LogNorm(vmin=Q_k.min(), vmax=Q_k.max()))

plt.tight_layout()
plt.show()

# Derivative regularization parameters
fig, axes = plt.subplots(1, 2, sharey='row', figsize=(8, 3.6))
fig.suptitle("Evolution of SOCC derivative regularization parameters")

Ap_k_image = axes[0].imshow(Ap_k, extent=(0, N_iter, w_points[-1], w_points[1]))
axes[0].set_title(r"$A'_k$")
axes[0].set_xlabel(r"# iteration")
axes[0].set_ylabel(r"$\omega_k$")
plt.colorbar(Ap_k_image, ax=axes[0], orientation='horizontal')

D_k_image = axes[1].imshow(D_k, extent=(0, N_iter, w_points[-1], w_points[1]))
axes[1].set_title(r"$D_k$")
axes[1].set_xlabel(r"# iteration")
axes[1].set_ylabel(r"$\omega_k$")
plt.colorbar(D_k_image,
             ax=axes[1],
             orientation='horizontal',
             norm=LogNorm(vmin=D_k.min(), vmax=D_k.max()))

plt.tight_layout()
plt.show()

# Second derivative regularization parameters
fig, axes = plt.subplots(1, 2, sharey='row', figsize=(8, 3.6))
fig.suptitle("Evolution of SOCC second derivative regularization parameters")

App_k_image = axes[0].imshow(App_k,
                             extent=(0, N_iter, w_points[-2], w_points[1]))
axes[0].set_title(r"$A''_k$")
axes[0].set_xlabel(r"# iteration")
axes[0].set_ylabel(r"$\omega_k$")
plt.colorbar(App_k_image, ax=axes[0], orientation='horizontal')

B_k_image = axes[1].imshow(B_k, extent=(0, N_iter, w_points[-2], w_points[1]))
axes[1].set_title(r"$D_k$")
axes[1].set_xlabel(r"# iteration")
axes[1].set_ylabel(r"$\omega_k$")
plt.colorbar(B_k_image,
             ax=axes[1],
             orientation='horizontal',
             norm=LogNorm(vmin=B_k.min(), vmax=B_k.max()))

plt.tight_layout()
plt.show()
