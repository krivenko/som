# This plotting script is largely based on a work of Malte Harland
# mharland@physnet.uni-hamburg.de

from pytriqs.gf.local import *
from pytriqs.gf.local.descriptors import *
from pytriqs.archive import HDFArchive
from pytriqs.statistics.histograms import Histogram
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

arch = HDFArchive('semicircle.h5','r')
pp = PdfPages('semicircle.pdf')

D = arch['D']
abs_errors = arch['abs_errors']

def make_ref(g):
    g << SemiCircular(D)

def plot_A_w(g_w, g_w_ref, fig):
    ax = fig.add_axes([.1,.55,.85,.4])
    w_mesh = [w for w in g_w.mesh]
    ax.plot(w_mesh, -g_w.data[:,0,0].imag/np.pi, color = 'red', linewidth = 0.6, label = 'SOM')
    ax.plot(w_mesh, -g_w_ref.data[:,0,0].imag/np.pi, color = 'blue', linewidth = 0.6, linestyle='dashed', label = 'reference')
    ax.legend()
    ax.set_xlabel("$\\omega$")
    ax.set_ylabel("$A(\\omega)$")
    ax.set_xlim(w_mesh[0].real, w_mesh[-1].real)
    ax.set_ylim(0, 0.3*np.pi/D)

def plot_histogram(histos, fig):
    hist = histos[0]
    ax = fig.add_axes([.3+2*.65/3.,.1,.65/3.,.35])
    dx = (hist.limits[1] - hist.limits[0]) / len(hist)
    x = np.linspace(hist.limits[0], hist.limits[1], len(hist.data))
    y = hist.data
    ax.bar(x, y, dx, color = 'green')
    ax.set_xlabel("$D$")
    ax.set_ylabel("$P(D)$")
    ax.set_xlim(*hist.limits)
    ax.set_ylim(bottom=0)

def make_g_iw_page(gr, s):
    g_w_ref = gr['g_w'].copy()
    make_ref(g_w_ref)

    # New figure
    fig = plt.figure()
    fig.suptitle("abs_error = %.4f, $G(i\\omega)$" % s)

    # Plot DOS
    plot_A_w(gr['g_w'], g_w_ref, fig)
    # Plot histogram
    plot_histogram(gr['histograms'], fig)
    # Plot Im G(i\omega) and Im G_{rec}(i\omega)
    g_iw = gr['g']
    g_iw_rec = gr['g_rec']
    iw_mesh = [w.imag for w in g_iw.mesh if w.imag > 0]
    ax = fig.add_axes([.1,.1,.54,.35])
    ax.plot(iw_mesh ,g_iw.data[len(iw_mesh):,0,0].imag, color = 'blue', linewidth = 0.6, label = 'original')
    ax.plot(iw_mesh ,g_iw_rec.data[len(iw_mesh):,0,0].imag, color = 'red', linewidth = 0.6, label = 'reconstructed')
    ax.set_xlabel("$\\omega_n$")
    ax.set_ylabel("$\\Im G(i\\omega)$")
    ax.legend()

    pp.savefig(fig)

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
    tau_mesh = [tau for tau in g_tau.mesh]
    ax = fig.add_axes([.1,.1,.54,.35])
    ax.plot(tau_mesh, g_tau.data[:,0,0], color = 'blue', linewidth = 0.6, label = 'original')
    ax.plot(tau_mesh, g_tau_rec.data[:,0,0], color = 'red', linewidth = 0.6, label = 'reconstructed')
    ax.set_xlim(0, g_tau.beta)
    ax.set_ylim(-1, 0)
    ax.set_xlabel("$\\tau$")
    ax.set_ylabel("$G(\\tau)$")
    ax.legend()

    pp.savefig(fig)

def make_g_l_page(gr, s):
    g_w_ref = gr['g_w'].copy()
    make_ref(g_w_ref)

    # New figure
    fig = plt.figure()
    fig.suptitle("abs_error = %.4f, $G(\\ell)$" % s)

    # Plot DOS
    plot_A_w(gr['g_w'], g_w_ref, fig)
    # Plot histogram
    plot_histogram(gr['histograms'], fig)
    # Plot G(\ell) and Im G_{rec}(\ell)
    g_l = gr['g']
    g_l_rec = gr['g_rec']
    l_mesh = [l for l in g_l.mesh]
    ax = fig.add_axes([.1,.1,.54,.35])
    ax.plot(l_mesh, g_l.data[:,0,0], color = 'blue', linewidth = 0.6, label = 'original')
    ax.plot(l_mesh, g_l_rec.data[:,0,0], color = 'red', linewidth = 0.6, label = 'reconstructed')
    ax.set_xlim(0, l_mesh[-1].real)
    ax.set_xlabel("$\\ell$")
    ax.set_ylabel("$G(\\ell)$")
    ax.legend()

    pp.savefig(fig)

for s in abs_errors:
    abs_err_gr = arch['abs_error_%.4f' % s]
    make_g_iw_page(abs_err_gr['g_iw'], s)
    make_g_tau_page(abs_err_gr['g_tau'], s)
    make_g_l_page(abs_err_gr['g_l'], s)

pp.close()
