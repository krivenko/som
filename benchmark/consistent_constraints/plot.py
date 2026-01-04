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

from h5 import HDFArchive
from triqs.gf import Gf                                             # noqa: F401
from triqs.gf.descriptors import SemiCircular
from triqs.plot.mpl_interface import oplot
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from sys import argv

arch_name = argv[1]

arch = HDFArchive(arch_name, 'r')
pp = PdfPages(arch_name.removesuffix(".h5") + '.pdf')

D = arch['input']['D']
energy_window = arch['cc_update_off']['accumulate_params']['energy_window']

fig = plt.figure()
fig.suptitle("Spectral functions")

# Plot referenced spectral function
g_w_ref = arch['cc_update_off']['som_output']['g_w'].copy()
g_w_ref << SemiCircular(D)
oplot(g_w_ref, mode='S', lw=1.0, label="Reference")

# Plot default model
w = np.array([float(w) for w in g_w_ref.mesh])
default_model = arch['cc_update_off']['socc_output']['params']['default_model']
plt.plot(w, default_model, lw=0.5, ls="--", label="SOCC default model")

for cc_update_gr_name in ('cc_update_off', 'cc_update_on'):
    gr = arch[cc_update_gr_name]

    g_w_som = gr['som_output']['g_w']
    chi_som = gr['som_output']['chi2'][0]

    g_w_socc = gr['socc_output']['g_w']
    chi_socc = gr['socc_output']['chi2'][0]

    cc_label = "CC+" if cc_update_gr_name == 'cc_update_on' else ''

    # Plot SOM spectral function
    oplot(g_w_som,
          mode='S',
          lw=0.5,
          label="%sSOM, $\\chi^2=%f$" % (cc_label, chi_som))

    # Plot SOCC spectral function
    oplot(g_w_socc,
          mode='S',
          lw=0.5,
          label="%sSOCC, $\\chi^2=%f$" % (cc_label, chi_socc))

plt.xlim(energy_window)
plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")
plt.legend(loc="upper left", prop={"size": 8})

pp.savefig(fig)

#
# Plot SOCC regularization parameters
#

for cc_update_gr_name in ('cc_update_off', 'cc_update_on'):
    cc_update_tile = "CC update " + \
        ("enabled" if cc_update_gr_name == "cc_update_on" else "disabled")

    gr = arch[cc_update_gr_name]

    socc_iterations = gr['socc_output']['iterations']
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
    fig, axes = plt.subplots(1, 2, sharey='row')
    fig.suptitle("Convergence of SOCC amplitude regularization parameters,\n%s"
                 % cc_update_tile)

    A_k_image = axes[0].imshow(A_k, extent=(0, N_iter, w[-1], w[0]))
    axes[0].set_title("$A_k$")
    axes[0].set_xlabel("# iteration")
    axes[0].set_ylabel(r"$\omega_k$")
    plt.colorbar(A_k_image, ax=axes[0], orientation='horizontal')

    Q_k_image = axes[1].imshow(Q_k, extent=(0, N_iter, w[-1], w[0]))
    axes[1].set_title("$Q_k$")
    axes[1].set_xlabel("# iteration")
    axes[1].set_ylabel(r"$\omega_k$")
    plt.colorbar(Q_k_image,
                 ax=axes[1],
                 orientation='horizontal',
                 norm=LogNorm(vmin=Q_k.min(), vmax=Q_k.max()))

    pp.savefig(fig)

    # Derivative regularization parameters

    fig, axes = plt.subplots(1, 2, sharey='row')
    fig.suptitle("Convergence of SOCC derivative regularization parameters,\n%s"
                 % cc_update_tile)

    Ap_k_image = axes[0].imshow(Ap_k, extent=(0, N_iter, w[-1], w[1]))
    axes[0].set_title("$A'_k$")
    axes[0].set_xlabel("# iteration")
    axes[0].set_ylabel(r"$\omega_k$")
    plt.colorbar(Ap_k_image, ax=axes[0], orientation='horizontal')

    D_k_image = axes[1].imshow(D_k, extent=(0, N_iter, w[-1], w[1]))
    axes[1].set_title("$D_k$")
    axes[1].set_xlabel("# iteration")
    axes[1].set_ylabel(r"$\omega_k$")
    plt.colorbar(D_k_image,
                 ax=axes[1],
                 orientation='horizontal',
                 norm=LogNorm(vmin=D_k.min(), vmax=D_k.max()))

    pp.savefig(fig)

    # Second derivative regularization parameters

    fig, axes = plt.subplots(1, 2, sharey='row')
    fig.suptitle(
        "Convergence of SOCC second derivative regularization parameters,\n%s"
        % cc_update_tile)

    App_k_image = axes[0].imshow(App_k, extent=(0, N_iter, w[-2], w[1]))
    axes[0].set_title("$A''_k$")
    axes[0].set_xlabel("# iteration")
    axes[0].set_ylabel(r"$\omega_k$")
    plt.colorbar(App_k_image, ax=axes[0], orientation='horizontal')

    B_k_image = axes[1].imshow(B_k, extent=(0, N_iter, w[-2], w[1]))
    axes[1].set_title("$D_k$")
    axes[1].set_xlabel("# iteration")
    axes[1].set_ylabel(r"$\omega_k$")
    plt.colorbar(B_k_image,
                 ax=axes[1],
                 orientation='horizontal',
                 norm=LogNorm(vmin=B_k.min(), vmax=B_k.max()))

    pp.savefig(fig)

pp.close()
