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

from h5 import HDFArchive
from triqs.gf import Gf                                             # noqa: F401
from matplotlib import pyplot as plt
from triqs.plot.mpl_interface import oplot
from matplotlib.backends.backend_pdf import PdfPages

arch = HDFArchive('bosonautocorr.h5', 'r')
pp = PdfPages('bosonautocorr.pdf')

chi_tau = arch["BosonCorr"]["input"]
chi_w = arch["BosonCorr"]["output"]
chi_w_symm = arch["BosonAutoCorr"]["output"]
chi_rec_tau = arch["BosonCorr"]["rec"]
chi_rec_tau_symm = arch["BosonAutoCorr"]["rec"]

tail = arch["BosonCorr"]["tail"]
tail_symm = arch["BosonAutoCorr"]["tail"]

oplot(chi_w[0, 0], mode='R', lw=0.5, label="BosonCorr, Re")
oplot(chi_w[0, 0], mode='I', lw=0.5, label="BosonCorr, Im")
oplot(chi_w_symm[0, 0], mode='R', lw=0.5, label="BosonAutoCorr, Re")
oplot(chi_w_symm[0, 0], mode='I', lw=0.5, label="BosonAutoCorr, Im")

w_mesh = list(map(float, chi_w.mesh))
plt.xlim((w_mesh[0], w_mesh[-1]))
plt.xlabel("$\\omega$")
plt.ylabel("$\\chi(\\omega)$")
plt.legend()
pp.savefig(plt.gcf())
plt.clf()

plt.plot(tail.flatten().real, lw=0.5, label="tail coefficients, BosonCorr")
plt.plot(tail_symm.flatten().real,
         lw=0.5,
         label="tail coefficients, BosonAutoCorr")
plt.xlim((0, len(tail) - 1))
plt.xlabel("order")
plt.legend()
pp.savefig(plt.gcf())
plt.clf()

oplot(chi_tau[0, 0], mode='R', lw=0.5, label="input")
oplot(chi_rec_tau[0, 0], mode='R', lw=0.5, label="rec, BosonCorr")
oplot(chi_rec_tau_symm[0, 0], mode='R', lw=0.5, label="rec, BosonAutoCorr")

tau_mesh = list(map(float, chi_tau.mesh))
plt.xlim((tau_mesh[0], tau_mesh[-1]))
plt.xlabel("$\\tau$")
plt.ylabel("$\\chi(\\tau)$")
plt.legend()
pp.savefig(plt.gcf())

pp.close()
