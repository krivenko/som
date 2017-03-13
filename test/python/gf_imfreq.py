##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016 by I. Krivenko
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

from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from pytriqs.applications.analytical_continuation.som import Som
from pytriqs.utility.comparison_tests import *

arch = HDFArchive('gf_imfreq.ref.h5', 'r')

g_iw = arch['g_iw']

cont = Som(g_iw, kind = "FermionGf")

run_params = {'energy_window' : (-5,5)}
run_params['verbosity'] = 3
run_params['random_seed'] = 34788 + 928374 * mpi.rank
run_params['adjust_f'] =  False
run_params['adjust_l'] = False
run_params['t'] = 100
run_params['f'] = 50
run_params['l'] = 30
run_params['make_histograms'] = True

cont.run(**run_params)

g_rec_iw = g_iw.copy()
g_rec_iw << cont

g_w = GfReFreq(window = (-6.0,6.0), n_points = 1200, indices = g_iw.indices)
g_w << cont

if mpi.is_master_node():
#    del arch
#    with HDFArchive('gf_imfreq.ref.h5', 'a') as arch:
#        arch['g_rec_iw'] = g_rec_iw
#        arch['g_w'] = g_w
#        arch['histograms'] = cont.histograms
    assert_gfs_are_close(g_rec_iw, arch['g_rec_iw'])
    assert_gfs_are_close(g_w, arch['g_w'], 6e-4)
    assert_arrays_are_close(g_w.tail.data, arch['g_w'].tail.data, 6e-4)
    assert_arrays_are_close(cont.histograms[0].data, arch['histograms'][0].data)
    assert_arrays_are_close(cont.histograms[1].data, arch['histograms'][1].data)
