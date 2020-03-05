##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#
# Compare results of this benchmark from SOM 1.x
#

from pytriqs.gf import *
from pytriqs.gf.descriptors import *
from pytriqs.archive import HDFArchive
from som import Som
from scipy.integrate import quad
import numpy as np
import h5py
from itertools import product

atol = 1e-12

arch_name = "all_kernels.np1.h5"
arch_old_name = "all_kernels.np1.1.x.h5"

arch = HDFArchive(arch_name, 'r')
arch_old = HDFArchive(arch_old_name, 'r')

for kind in ("FermionGf", "BosonCorr", "BosonAutoCorr", "ZeroTemp"):
    kind_msg = kind + " kernels"
    print("=" * len(kind_msg))
    print(kind_msg)
    print("=" * len(kind_msg))

    if not kind in arch:
        print("WARNING: /%s group is not in %s, skipping..." %(kind, arch_name))
        continue
    if not kind in arch_old:
        print("WARNING: /%s group is not in %s, skipping..." %(kind, arch_old_name))
        continue

    for mesh in ("imtime", "imfreq", "legendre"):
        print("-"*len(mesh))
        print(mesh)
        print("-"*len(mesh))

        if not mesh in arch[kind]:
            print("WARNING: /%s/%s group is not in %s, skipping..." %
                  (kind, mesh, arch_name))
            continue
        if not mesh in arch_old[kind]:
            print("WARNING: /%s/%s group is not in %s, skipping..." %
                  (kind, mesh, arch_old_name))
            continue

        gr, gr_old = arch[kind][mesh], arch_old[kind][mesh]

        # Input
        inp, inp_old = gr["input"], gr_old["input"]
        assert inp.mesh == inp_old.mesh
        if np.allclose(inp.data, inp_old.data, atol=atol, rtol=0):
            print("Input matches")
        else:
            print("Input max deviation: %f" %
                np.max(np.abs(inp.data - inp_old.data)))

        # Solutions
        sols, sols_old = gr["solutions"], gr_old["solutions"]
        for n, (s, s_old) in enumerate(zip(sols, sols_old)):
            if len(s) != len(s_old):
                print("Solution %d: mismatching number of rectangles (%d vs %d)" %
                      (n, len(s), len(s_old)))
            else:
                if np.allclose(s, s_old, atol=atol, rtol=0):
                    print("Solution %d matches" % n)
                else:
                    for nr, (r, r_old) in enumerate(zip(s, s_old)):
                        if np.allclose(r, r_old, atol=atol, rtol=0):
                            continue
                        else:
                            print("Solution %d, rectangle %d (deviation): %s" %
                                  (n, nr, np.asarray(r) - np.asarray(r_old)))

        # Output
        out, out_old = gr["output"], gr_old["output"]
        assert out.mesh == out_old.mesh
        if np.allclose(out.data, out_old.data, atol=atol, rtol=0):
            print("Output matches")
        else:
            print("Output max deviation: %f" %
                np.max(np.abs(out.data - out_old.data)))

        # Reconstructed
        rec, rec_old = gr["rec"], gr_old["rec"]
        assert rec.mesh == rec_old.mesh
        if np.allclose(rec.data, rec_old.data, atol=atol, rtol=0):
            print("Reconstruction matches")
        else:
            print("Reconstruction max deviation: %f" %
                np.max(np.abs(rec.data - rec_old.data)))

        # Tail
        tail_mask = gr_old["output/singularity/mask"]
        tail = gr["output_tail"]
        tail_old = gr_old["output/singularity/data"][1:,:,:]
        for i, j in product(*map(range, tail_mask.shape)):
            m = tail_mask[i, j]+1
            if np.allclose(tail[:m, i, j], tail_old[:, i, j], atol=atol, rtol=0):
                print("Tail element (%d,%d) matches" % (i,j))
            else:
                print("Tail element (%d,%d) max deviation: %f" %
                      (i,j,np.max(np.abs(tail[:m, i, j] - tail_old[:, i, j]))))

        # Histograms
        hists, hists_old = gr["histograms"], gr_old["histograms"]
        for n, (h, h_old) in enumerate(zip(hists, hists_old)):
            if np.allclose(h.data, h_old.data, atol=atol, rtol=0):
                print("Histogram %d matches" % n)
            else:
                print("Histogram %d max deviation: %f" % \
                      (n, np.max(np.abs(h.data - h_old.data))))

