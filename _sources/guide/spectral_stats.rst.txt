.. _spectral_stats:

Statistical analysis of ensembles of spectral functions
=======================================================

.. currentmodule:: som.spectral_stats

A theoretical framework presented in Sections I-II of [GMPPS2017]_ offers a way
to estimate statistical properties (averages, dispersions and two-point
correlations) of ensembles of
:ref:`particular solutions <particular_solutions>`. It is implemented as a
handful of functions in the :mod:`som.spectral_stats` module.

The analysis is built around the notion of spectral integrals. Given a set of
finite intervals :math:`\{\Delta_m\}` centered around points :math:`z_m`,
the spectral integrals are defined as

.. math::

    i_m = \int_{-\infty}^\infty dz \bar K(m, z) A(z),

where kernel :math:`\bar K(m, z)` is a "resolution function" localized inside
:math:`\{\Delta_m\}`.

Function :func:`som.spectral_stats.spectral_integral` evaluates :math:`i_m` for
a single interval :math:`\Delta_m`, for all intervals induced by a
:class:`real frequency mesh <triqs.gf.meshes.MeshReFreq>`, or for an
arbitrary set of intervals. It supports a few kernels :math:`\bar K(m, z)`
switched by the keyword argument ``resolution_function``,

.. _resolution_functions:

- ``resolution_function=rectangle``:

    .. math::
        \bar K(m, z) = \frac{1}{|\Delta_m|}
            \theta\left(\frac{|\Delta_m|}{2} - |z - z_m|\right).

- ``resolution_function=lorentzian``:

    .. math::
        \bar K(m, z) = \frac{1}{\pi}
            \frac{|\Delta_m|/2}{(z - z_m)^2 + (|\Delta_m|/2)^2}.


- ``resolution_function=gaussian``:

    .. math::
        \bar K(m, z) = \frac{1}{\sqrt{2\pi(|\Delta_m|/2)^2}}
            \exp\left(-\frac{(z - z_m)^2}{2(|\Delta_m|/2)^2}\right).

Having accumulated a set (ensemble) of :math:`J` particular solutions, one can
study statistical properties of spectral integrals :math:`i_m^{(j)}` derived
from them.

Function :func:`som.spectral_stats.spectral_avg` computes averages over the
ensemble,

.. math::

    \bar{i}_m = \frac{1}{J} \sum_{j=1}^J i_m^{(j)}.

These averages can be interpreted as an estimate of the final spectrum.

While horizontal "error bars" on :math:`\bar{i}_m` are simply
:math:`|\Delta_m|/2`, the vertical ones (:math:`\sigma_m`) are estimated from
dispersions

.. math::

    \sigma^2_m = \frac{1}{J} \sum_{j=1}^J (i_m^{(j)} - \bar{i}_m)^2.

Spectral dispersions are computed by :func:`som.spectral_stats.spectral_disp`,
and :func:`som.spectral_stats.spectral_corr` returns a full matrix of
two-points correlators

.. math::

    \sigma_{mm'} = \frac{1}{J} \sum_{j=1}^J (i_m^{(j)} - \bar{i}_m)
                                            (i_{m'}^{(j)} - \bar{i}_{m'}).
