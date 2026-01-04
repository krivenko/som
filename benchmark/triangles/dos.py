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

import numpy as np
from scipy.integrate import quad
from triqs.gf.descriptors import Function


# Integral kernel
def make_kernel(beta):
    def K(tau, e):
        if e > 0:
            return -np.exp(-e*tau.real) / (1 + np.exp(-e*beta))
        else:
            return -np.exp(e*(beta-tau.real)) / (1 + np.exp(e*beta))
    return K


e_min = -0.9
e_low = 2.0
e_max = 6.0


# DOS
def A(w):
    if w >= -0.9 and w < 1.0:
        return (0.25*w + 0.225)/1.18875000078
    if w >= 1.0 and w < 2.0:
        return (-0.475*w + 0.95)/1.18875000078
    if w >= 4.0 and w < 4.5:
        return (w - 4.0)/1.18875000078
    if w >= 4.5 and w < 6.0:
        return (-1.0/3.0*w + 2.0)/1.18875000078
    return 0


def make_g_tau(g_tau):
    K = make_kernel(g_tau.mesh.beta)
    g_tau << Function(
        lambda tau: quad(lambda w: K(tau, w) * A(w), e_min, e_max, limit=100)[0]
    )
