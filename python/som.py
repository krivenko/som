##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2017 by I. Krivenko
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

from core import SomCore
import numpy as np

class Som(SomCore):
    """Stochastic Optimization Method"""

    def __init__(self, g, s = None, kind = "FermionGf", norms = np.array([])):
        if s is None:
            s = g.copy()
            s.data[:,Ellipsis] = np.eye(s.target_shape[0])
        if isinstance(norms,float) or isinstance(norms,int):
            norms = norms * np.ones((g.target_shape[0],))
        SomCore.__init__(self, g, s, kind, norms)
