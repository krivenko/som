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


# Theta-function
def Theta(x): return 0.5 * (np.sign(x) + 1)


# Integral kernel
def make_kernel(beta):
    return lambda tau, e: -np.exp(-e*tau.real) / (1 + np.exp(-e*beta))


# Delta-function
e_delta = 0.03      # position
Z_delta = 0.07      # weight
# Continuous DOS
e_th = 0.04             # microgap
e_gap = e_th - e_delta
e_con = 0.566

e_min = 0.0
e_max = 2.0


def A_delta(w):
    d = 1e-5
    return Z_delta * d / ((w - e_delta)**2 + d**2)


def A_con(w):
    if w >= e_th and w <= e_con:
        return Z_delta * np.sqrt(w - e_th) / \
               (2*np.pi * np.sqrt(e_gap)*((w - e_th) + e_gap))
    else:
        return 0


def A_triangle(w):
    res = 0
    if w >= e_con and w < 1:
        f = A_con(e_con)
        res += (2.30415 - 2.30415*f)*w + (-1.30415 + 2.30415*f)
    elif w >= 1.0 and w < 2.0:
        res += (2.0 - w)
    return 0.868*res


def make_g_tau(g_tau):
    K = make_kernel(g_tau.mesh.beta)
    g_tau << Function(lambda tau: Z_delta*K(tau, e_delta) +
                      quad(lambda w: K(tau, w)*(A_con(w) + A_triangle(w)),
                           e_min, e_max)[0])


def A(w): return A_delta(w) + A_con(w) + A_triangle(w)

# if   w>=0.5 and w<1:       res += 0.533615 * (4.0*w-2.0)
# elif w>=1.0 and w<2.0:     res += 0.533615 * (4.0-2.0*w)
