##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
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

from typing import List

import cmath
import numpy as np

from triqs.gf import (Gf,
                      GfImFreq,
                      MeshImFreq,
                      MeshImTime,
                      MeshLegendre,
                      Fourier)
from triqs.stat import Histogram

from .som_core import SomCore


class Som(SomCore):
    """Implementation of the Stochastic Optimization Method."""

    def __init__(self,
                 rhs: Gf,
                 errors: Gf,
                 kind: str = "FermionGf",
                 norms=None,
                 *,
                 filtering_levels=None,
                 ):
        r"""
        :param rhs: Right hand side of the :ref:`integral equation
                  <integral_equation>` to be solved, defined on
                  :class:`triqs.gf.meshes.MeshImTime`,
                  :class:`triqs.gf.meshes.MeshImFreq` or
                  :class:`triqs.gf.meshes.MeshLegendre`.
                  The target shape of ``rhs`` must be :math:`M{\times}M`.
                  If :math:`M>1`, only its diagonal matrix elements will be
                  considered and used as input data for :math:`M` independent
                  continuation problems.
        :type rhs: :class:`triqs.gf.gf.Gf`

        :param errors: Either :ref:`error bars <error_bars>` :math:`\sigma_n`
                       (GF container of the same type and target shape as
                       ``rhs``) or a packed list of :ref:`covariance matrices
                       <cov_matrix>` (``Gf(mesh=MeshProduct(rhs.mesh, rhs.mesh),
                       target_shape=[M])`` with each target element
                       corresponding to one covariance matrix).

        :param kind: Selection of the :ref:`observable kind <observables>`
                     (and its respective integral kernel) to be used. One of
                     ``FermionGf``, ``FermionGfSymm``, ``BosonCorr``,
                     ``BosonAutoCorr``, ``ZeroTemp``. Defaults to ``FermionGf``.
        :type kind: :class:`str`, optional

        :param norms: Requested :ref:`solution norms <solution_norm>` either as
                      a list of :math:`M` real numbers or as a single real
                      number for all :math:`M` continuation problems. For
                      observable kinds ``FermionGf``, ``FermionGfSymm`` and
                      ``ZeroTemp`` this parameter is optional and defaults to
                      1.0. For ``BosonCorr`` and ``BosonAutoCorr`` it must be
                      provided by the user.

                      .. seealso::

                        When unknown a priori in the latter case, the norms can
                        be approximately estimated from ``rhs`` using
                        :func:`estimate_boson_corr_spectrum_norms`.

        :type norms: :class:`float` or :class:`list` [:class:`float`]

        :param filtering_levels: :ref:`Filtering levels <cov_matrix_filtered>`
                                 for covariance matrices either as
                                 a list of :math:`M` real numbers or as a
                                 single real number for all :math:`M`
                                 continuation problems. Can only be specified
                                 when the covariance matrices are used.
                                 Defaults to 0.
        :type filtering_levels: :class:`float` or
                                :class:`list` [:class:`float`], optional
        """

        if rhs.target_rank != 2:
            raise RuntimeError("'rhs' must be a square matrix-valued GF object")

        if norms is None:
            try:
                norms_ = {"FermionGf": 1.0,
                          "FermionGfSymm": 1.0,
                          "ZeroTemp":  1.0}[kind] * np.ones(rhs.target_shape[0])
            except KeyError:
                raise RuntimeError("A list of solution norms must be provided "
                                   "for observable kind " + kind)
        elif isinstance(norms, float) or isinstance(norms, int):
            norms_ = norms * np.ones(rhs.target_shape[0])
        else:
            norms_ = np.array(norms)

        # First, try to construct with the covariance matrices
        if isinstance(errors, Gf) \
                and errors.rank == 2 and errors.target_rank == 1:

            if filtering_levels is None:
                fl = np.array([])
            elif isinstance(filtering_levels, float) or \
                    isinstance(filtering_levels, int):
                fl = filtering_levels * np.ones(rhs.target_shape[0])
            else:
                fl = np.array(filtering_levels)

            SomCore.__init__(self, rhs, errors, kind, norms_, fl)

        # Then try the error bars
        elif isinstance(errors, Gf) \
                and errors.rank == 1 and errors.target_rank == 2:

            if filtering_levels is not None:
                raise RuntimeError(
                    "Argument 'filtering_levels' is accepted only when full "
                    "covariance matrices are provided")

            SomCore.__init__(self, rhs, errors, kind, norms_)

        # Give up
        else:
            raise RuntimeError("Argument 'errors' has unsupported format")


def estimate_boson_corr_spectrum_norms(chi: Gf) -> List[float]:
    r"""
    Given a correlator :math:`\chi` of bosonic or boson-like
    (fermion-number-conserving) operators, estimates the corresponding spectrum
    normalization constants :math:`\mathcal{N}`.

    Depending on the mesh :math:`\chi` is defined on, one of the following
    expressions is used.

    - Imaginary frequencies: :math:`\mathcal{N} = \pi\chi(i\Omega=0)`.
    - Imaginary time: :math:`\mathcal{N} = \pi\int_0^\beta d\tau \chi(\tau)`.
    - Legendre polynomial basis coefficients: :math:`\mathcal{N} =
      \pi\chi(\ell=0)`.

    :param chi: The correlator :math:`\chi`.
    :type chi: :class:`triqs.gf.gf.Gf`
    :return: A list of estimated spectrum normalization constants,
             one constant per diagonal matrix element of :math:`\chi`.
    :rtype: :class:`list` [:class:`float`]
    """
    assert isinstance(chi, Gf), "Expected a Green's function object"

    if chi.mesh.statistic != "Boson":
        raise ValueError("Wrong mesh statistics for bosonic correlator 'chi'")

    N = chi.target_shape[0]

    if isinstance(chi.mesh, MeshImFreq):
        W0 = list(chi.mesh.values()).index(0j)
        return np.pi * np.array([chi.data[W0, n, n].real for n in range(N)])

    elif isinstance(chi.mesh, MeshImTime):
        chi_iw = GfImFreq(beta=chi.mesh.beta,
                          statistic="Boson",
                          n_points=1,  # We need only the zero frequency
                          target_shape=chi.target_shape)
        chi_iw << Fourier(chi)
        return np.pi * np.array([chi_iw.data[0, n, n].real for n in range(N)])

    elif isinstance(chi.mesh, MeshLegendre):
        # \chi(i\Omega = 0) = \chi(l = 0)
        return np.pi * np.array([chi.data[0, n, n].real for n in range(N)])
    else:
        raise TypeError("Unexpected type of 'chi'")


def count_good_solutions(hist: Histogram,
                         good_chi_rel: float = 2.0,
                         good_chi_abs: float = cmath.inf) -> int:
    r"""
    Given a histogram of values :math:`\chi` for the objective function
    :math:`\chi^2`, counts the number of such solutions that
    :math:`\chi \leq` ``good_chi_abs`` and
    :math:`\chi/\chi_\mathrm{min} \leq` ``good_chi_rel``.

    :param hist: Histogram to analyze.
    :type hist: :class:`triqs.stat.histograms.Histogram`
    :param good_chi_rel: :math:`\chi/\chi_\mathrm{min}` threshold for good
                         solutions.
    :type good_chi_rel: :class:`float`, optional
    :param good_chi_abs: :math:`\chi` threshold for good solutions.
    :type good_chi_abs: :class:`float`, optional
    :rtype: :class:`int`
    """
    chi_min = hist.limits[0]
    chi_c = min(good_chi_abs, chi_min * good_chi_rel)
    return int(sum(c for n, c in enumerate(hist.data)
                   if hist.mesh_point(n) <= chi_c))
