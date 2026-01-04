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
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot
from matplotlib.backends.backend_pdf import PdfPages

arch = HDFArchive('fermiongfsymm.h5', 'r')
pp = PdfPages('fermiongfsymm.pdf')

g_tau = arch["FermionGf"]["input"]
g_w = arch["FermionGf"]["output"]
g_w_symm = arch["FermionGfSymm"]["output"]
g_rec_tau = arch["FermionGf"]["rec"]
g_rec_tau_symm = arch["FermionGfSymm"]["rec"]

tail = arch["FermionGf"]["tail"]
tail_symm = arch["FermionGfSymm"]["tail"]

oplot(g_w[0, 0], mode='R', lw=0.5, label="FermionGf, Re")
oplot(g_w[0, 0], mode='I', lw=0.5, label="FermionGf, Im")
oplot(g_w_symm[0, 0], mode='R', lw=0.5, label="FermionGfSymm, Re")
oplot(g_w_symm[0, 0], mode='I', lw=0.5, label="FermionGfSymm, Im")

w_mesh = list(map(float, g_w.mesh))
plt.xlim((w_mesh[0], w_mesh[-1]))
plt.xlabel("$\\omega$")
plt.ylabel("$G(\\omega)$")
plt.legend()
pp.savefig(plt.gcf())
plt.clf()

plt.plot(tail.flatten().real, lw=0.5, label="tail coefficients, FermionGf")
plt.plot(tail_symm.flatten().real, lw=0.5,
         label="tail coefficients, FermionGfSymm")
plt.xlim((0, len(tail) - 1))
plt.xlabel("order")
plt.legend()
pp.savefig(plt.gcf())
plt.clf()

oplot(g_tau[0, 0], mode='R', lw=0.5, label="input")
oplot(g_rec_tau[0, 0], mode='R', lw=0.5, label="rec, FermionGf")
oplot(g_rec_tau_symm[0, 0], mode='R', lw=0.5, label="rec, FermionGfSymm")

tau_mesh = list(map(float, g_tau.mesh))
plt.xlim((tau_mesh[0], tau_mesh[-1]))
plt.xlabel("$\\tau$")
plt.ylabel("$G(\\tau)$")
plt.legend()
pp.savefig(plt.gcf())

pp.close()
