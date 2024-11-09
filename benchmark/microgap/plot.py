##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2024 Igor Krivenko
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

# This plotting script is largely based on a work of Malte Harland
# mharland@physnet.uni-hamburg.de

from h5 import HDFArchive
from triqs.gf import Gf                                             # noqa: F401
from triqs.stat.histograms import Histogram                         # noqa: F401
from triqs.gf.descriptors import Function
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.integrate import quad

from dos import A, e_min, e_max, e_con

arch = HDFArchive('microgap.h5', 'r')
pp = PdfPages('microgap.pdf')

abs_errors = arch['abs_errors']


def make_ref(g):
    print("Norm of the reference DOS:",
          quad(A, e_min, e_max, limit=100)[0].real)
    g << Function(lambda w: -1j*np.pi * A(w.real))


def plot_A_w(g_w, g_w_ref, fig):
    w_mesh = [w.real for w in g_w.mesh]

    ax = fig.add_axes([.1, .54, .55, .4])
    ax.plot(w_mesh,
            -g_w.data[:, 0, 0].imag/np.pi,
            color='red',
            linewidth=0.6,
            label='SOM')
    ax.plot(w_mesh,
            -g_w_ref.data[:, 0, 0].imag/np.pi,
            color='blue',
            linewidth=0.6,
            linestyle='dashed',
            label='reference')
    ax.legend(fontsize=10.0)
    ax.set_xlabel("$\\omega$")
    ax.set_ylabel("$A(\\omega)$")
    ax.set_xlim(e_min, e_max)
    ax.set_ylim(0, 1.2)

    ax = fig.add_axes([.3+2*.65/3., .55, .65/3, .4])
    ax.plot(w_mesh,
            -g_w.data[:, 0, 0].imag/np.pi,
            color='red',
            linewidth=0.6,
            label='SOM')
    ax.plot(w_mesh,
            -g_w_ref.data[:, 0, 0].imag/np.pi,
            color='blue',
            linewidth=0.6,
            linestyle='dashed',
            label='reference')
    ax.legend(fontsize=10.0)
    ax.set_xlabel("$\\omega$")
    ax.set_ylabel("$A(\\omega)$")
    ax.set_xlim(e_min, e_con)
    ax.set_ylim(0, 1.2)


def plot_histogram(histos, fig):
    hist = histos[0]
    ax = fig.add_axes([.3+2*.65/3., .1, .65/3., .35])
    dx = (hist.limits[1] - hist.limits[0]) / len(hist)
    x = np.linspace(hist.limits[0], hist.limits[1], len(hist.data))
    y = hist.data
    ax.bar(x, y, dx, color='green')
    ax.set_xlabel(r"$\chi$")
    ax.set_ylabel(r"$P(\chi)$")
    ax.set_xlim(*hist.limits)
    ax.set_ylim(bottom=0)


def make_g_tau_page(gr, s):
    g_w_ref = gr['g_w'].copy()
    make_ref(g_w_ref)

    # New figure
    fig = plt.figure()
    fig.suptitle("abs_error = %.4f, $G(\\tau)$" % s)

    # Plot DOS
    plot_A_w(gr['g_w'], g_w_ref, fig)
    # Plot histogram
    plot_histogram(gr['histograms'], fig)
    # Plot G(i\tau) and Im G_{rec}(i\tau)
    g_tau = gr['g']
    g_tau_rec = gr['g_rec']
    tau_mesh = [tau.real for tau in g_tau.mesh]
    ax = fig.add_axes([.1, .1, .54, .35])
    ax.plot(tau_mesh,
            g_tau.data[:, 0, 0].real,
            color='blue',
            linewidth=0.6,
            label='original')
    ax.plot(tau_mesh,
            g_tau_rec.data[:, 0, 0].real,
            color='red',
            linewidth=0.6,
            label='reconstructed')
    ax.set_xlim(0, g_tau.mesh.beta)
    ax.set_ylim(-1, 0)
    ax.set_xlabel("$\\tau$")
    ax.set_ylabel("$G(\\tau)$")
    ax.legend()

    pp.savefig(fig)


for s in abs_errors:
    abs_err_gr = arch['abs_error_%.4f' % s]
    make_g_tau_page(abs_err_gr['g_tau'], s)

pp.close()
