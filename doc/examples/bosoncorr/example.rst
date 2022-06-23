.. _example_bosoncorr:

Green’s function of bosons or transverse magnetic susceptibility
================================================================

A general correlation function of boson-like operators :math:`\hat O` and
:math:`\hat O^\dagger` is defined as

    .. math::

        \chi(\tau) =
            \langle\mathbb{T}_\tau \hat O(\tau) \hat O^\dagger(0)\rangle.

:math:`\chi(\tau)` is subject to the periodicity condition

    .. math::

        \chi(\tau+\beta) = \chi(\tau).

Examples of such correlators are the transverse magnetic susceptibility
:math:`\langle\mathbb{T}_\tau \hat S_+(\tau) \hat S_-(0)\rangle` and the
Green's function of bosons
:math:`\langle\mathbb{T}_\tau \hat b(\tau) \hat b^\dagger(0)\rangle`.

The auxiliary spectral function :math:`A(\epsilon) = \Im\chi(\epsilon)/\epsilon`
is non-negative but not necessarily symmetric.

.. warning::

    When :math:`\hat O = \hat O^\dagger` (for instance, in the case of charge or
    longitudinal magnetic susceptibility), it is strongly recommended to use the
    :ref:`BosonAutoCorr <bosonautocorr>` observable kind instead
    of :ref:`BosonCorr <bosoncorr>` shown here.

Perform analytic continuation using the :ref:`BosonCorr <bosoncorr>` kernel
---------------------------------------------------------------------------

.. literalinclude:: example.py

Download input file :download:`example.h5`.

Plot input and reconstructed correlators at Matsubara frequencies
-----------------------------------------------------------------

.. plot:: examples/bosoncorr/plot_chi_iw.py
    :include-source:
    :scale: 100

Plot the correlator on the real frequency axis and its tail coefficients
------------------------------------------------------------------------

.. plot:: examples/bosoncorr/plot_chi_w.py
    :include-source:
    :scale: 100

Plot :math:`\chi`-histograms
----------------------------

.. plot:: examples/bosoncorr/plot_hist.py
    :include-source:
    :scale: 100
