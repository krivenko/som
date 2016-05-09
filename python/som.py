################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2016 by I. Krivenko
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

r"""
DOC

"""
from core import SomCore
import numpy as np

class Som(SomCore):
    """Stochastic Optimization Method"""

    def __init__(self, g, s = None, kind = "FermionGf", norm = 1.0):
        if s is None:
            s = g.copy()
            s.data[:,Ellipsis] = np.eye(s.target_shape[0])
        SomCore.__init__(self, g, s, kind, norm)

    def __rlshift__(self, g):
        SomCore.__call__(self, g)
