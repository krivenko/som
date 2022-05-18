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
from triqs.gf import Gf
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

def count_good_solutions(hist, upper_lim=1):
    r"""
    Given a histogram of \chi-values,
    count the number of solutions with \chi/\chi_{min} <= 1 + upper_lim
    """
    d_max = hist.limits[0] * (1 + upper_lim)
    return int(sum(c for n, c in enumerate(hist.data)
                   if hist.mesh_point(n) <= d_max))
