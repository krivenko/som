##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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

from h5 import *
from triqs.gf import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from triqs.plot.mpl_interface import oplot

arch = HDFArchive('binning.h5', 'r')
pp = PdfPages('binning.pdf')

def plot_binning_nobinning(kind, ylabel):
    g_w = arch[kind]["output"]
    g_w_nobinning = arch[kind]["output_nobinning"]

    n_w_nobinning = len(g_w_nobinning.mesh)
    n_w = len(g_w.mesh)

    plt.clf()
    oplot(g_w_nobinning[0, 0],
          name = "without binning, $N_{\\omega} = %d$" % n_w_nobinning)
    oplot(g_w[0, 0], name = "with binning, $N_{\\omega} = %d$" % n_w)

    w_mesh = list(map(float, g_w.mesh))

    plt.title(f"{kind} kernel")
    plt.xlim((w_mesh[0], w_mesh[-1]))
    plt.xlabel("$\\omega$")
    plt.ylabel(ylabel)
    plt.legend()

    pp.savefig(plt.gcf())

plot_binning_nobinning("FermionGf", "$G(\\omega)$")
plot_binning_nobinning("BosonCorr", "$\\chi(\\omega)$")
plot_binning_nobinning("BosonAutoCorr", "$\\chi(\\omega)$")
plot_binning_nobinning("ZeroTemp", "$G(\\omega)$")

pp.close()
