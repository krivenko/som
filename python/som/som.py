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
"""
Main module of SOM
"""

from .som_core import SomCore
from triqs.gf import (Gf,
                      GfImFreq,
                      GfImTime,
                      GfLegendre,
                      Fourier,
                      LegendreToMatsubara)
import numpy as np


class Som(SomCore):
    """Stochastic Optimization Method"""

    def __init__(self,
                 g,
                 errors,
                 kind="FermionGf",
                 norms=None,
                 *,
                 filtration_levels=None
                 ):

        if norms is None:
            try:
                norms_ = {"FermionGf": 1.0,
                          "FermionGfSymm": 1.0,
                          "ZeroTemp":  1.0}[kind] * np.ones(g.target_shape[0])
            except KeyError:
                raise RuntimeError("A list of solution norms must be provided "
                                   "for observable kind " + kind)
        elif isinstance(norms, float) or isinstance(norms, int):
            norms_ = norms * np.ones(g.target_shape[0])
        else:
            norms_ = np.array(norms)

        # First, try to construct with the covariance matrices
        if isinstance(errors, Gf) \
                and errors.rank == 2 and errors.target_rank == 1:

            if filtration_levels is None:
                fl = np.array([])
            elif isinstance(filtration_levels, float) or \
                    isinstance(filtration_levels, int):
                fl = filtration_levels * np.ones(g.target_shape[0])
            else:
                fl = np.array(filtration_levels)

            SomCore.__init__(self, g, errors, kind, norms_, fl)

        # Then try the error bars
        elif isinstance(errors, Gf) \
                and errors.rank == 1 and errors.target_rank == 2:

            if filtration_levels is not None:
                raise RuntimeError(
                    "Argument 'filtration_levels' is accepted only when full "
                    "covariance matrices are provided")

            SomCore.__init__(self, g, errors, kind, norms_)

        # Give up
        else:
            raise RuntimeError("Argument 'errors' has unsupported format")


def extract_boson_corr_spectrum_norms(chi):
    """
    Given a correlator of boson-like operators :math:`\chi` defined on any
    supported mesh, returns a list of spectrum normalization constants
    :math:`\mathcal{N} = \pi \chi(i\Omega = 0)`.
    """
    assert isinstance(chi, Gf), "Expected a Green's function object"

    if chi.mesh.statistic != "Boson":
        raise ValueError("Wrong mesh statistics for bosonic correlator 'chi'")

    N = chi.target_shape[0]

    if isinstance(chi, GfImFreq):
        W0 = list(chi.mesh.values()).index(0j)
        return np.pi * np.array([chi.data[W0, n, n].real for n in range(N)])

    elif isinstance(chi, GfImTime):
        chi_iw = GfImFreq(beta=chi.mesh.beta,
                          statistic="Boson",
                          n_points=1, # We need only the zero frequency
                          target_shape=chi.target_shape)
        chi_iw << Fourier(chi)
        return np.pi * np.array([chi_iw.data[0, n, n].real for n in range(N)])

    elif isinstance(chi, GfLegendre):
        # \chi(i\Omega = 0) = \chi(l = 0)
        return np.pi * np.array([chi.data[0, n, n].real for n in range(N)])
    else:
        raise TypeError("Unexpected type of 'chi'")


def count_good_solutions(hist, good_chi_rel=2.0, good_chi_abs=np.inf):
    r"""
    Given a histogram of values :math:`\chi` for the objective function
    :math:`\chi^2`, count the number of solutions such that
    :math:`\chi \leq \mathtt{good_chi_abs}` and
    :math:`\chi/\chi_\mathrm{min} \leq \mathtt{good_chi_rel}`.
    """
    chi_min = hist.limits[0]
    chi_c = min(good_chi_abs, chi_min * good_chi_rel)
    return int(sum(c for n, c in enumerate(hist.data)
                   if hist.mesh_point(n) <= chi_c))
