
Example: Correlator of boson-like operators
===========================================

Observable kind: `BosonCorr`.

A general correlation function of boson-like operators :math:`\hat O` and :math:`\hat O^\dagger`
is defined as

    .. math::

        \chi(\tau) = \langle\mathbb{T}_\tau \hat O(\tau) \hat O^\dagger(0)\rangle.

:math:`\chi(\tau)` is subject to an additional periodicity condition

    .. math::

        \chi(\tau+\beta) = \chi(\tau).

Examples of such correlators are the transverse magnetic susceptibility
:math:`\langle\mathbb{T}_\tau \hat S_+(\tau) \hat S_-(0)\rangle` and the
Green's function of bosons :math:`\langle\mathbb{T}_\tau \hat b(\tau) \hat b^\dagger(0)\rangle`.

The auxiliary spectral function :math:`A(\epsilon) = \Im\chi(\epsilon)/\epsilon`
is non-negative but not necessarily symmetric.

.. warning::

    When :math:`\hat A = \hat B` (for instance, in case of charge or longitudinal magnetic
    susceptibility), it is strongly recommended to use `BosonAutoCorr` observable kind instead.

Run analytical continuation
---------------------------

.. literalinclude:: example.py

Download input file :download:`example.h5`.

Plot input and reconstructed correlators at Matsubara frequencies
-----------------------------------------------------------------

.. plot:: examples/bosoncorr/plot_chi_iw.py
    :include-source:
    :scale: 100

Plot the correlator on the real frequency axis
----------------------------------------------

.. plot:: examples/bosoncorr/plot_chi_w.py
    :include-source:
    :scale: 100
