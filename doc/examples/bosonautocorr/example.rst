.. _example_bosonautocorr:

Charge susceptibility, longitudinal magnetic susceptibility and optical conductivity
====================================================================================

Correlator of a Hermitian operator :math:`\hat O` with itself is defined as

    .. math::

        \chi(\tau) = \langle\mathbb{T}_\tau \hat O(\tau) \hat O(0)\rangle.

Common examples of such correlators are the longitudinal magnetic susceptibility
:math:`\langle\mathbb{T}_\tau \hat S_z(\tau) \hat S_z(0)\rangle`, the
charge susceptibility
:math:`\langle\mathbb{T}_\tau \hat N(\tau) \hat N(0)\rangle` and the optical
conductivity :math:`\langle \mathbb{T}_\tau \hat j(\tau) \hat j(0)\rangle`.

For the Hermitian operators, the auxiliary spectral function
:math:`A(\epsilon) = \Im\chi(\epsilon)/\epsilon` is non-negative and symmetric.

This kind of correlators is treated by the :ref:`BosonAutoCorr <bosonautocorr>`
kernels, which are faster and more robust than :ref:`BosonCorr <bosoncorr>`.

Perform analytic continuation using the :ref:`BosonAutoCorr <bosonautocorr>` kernel
-----------------------------------------------------------------------------------

.. literalinclude:: example.py

Download input file :download:`input.h5`.

Plot input and reconstructed correlators at Matsubara frequencies
-----------------------------------------------------------------

.. plot:: examples/bosonautocorr/plot_chi_iw.py
    :include-source:
    :scale: 100

Plot the correlator on the real frequency axis and its tail coefficients
------------------------------------------------------------------------

.. plot:: examples/bosonautocorr/plot_chi_w.py
    :include-source:
    :scale: 100


Plot :math:`\chi`-histograms
----------------------------

.. plot:: examples/bosonautocorr/plot_hist.py
    :include-source:
    :scale: 100
