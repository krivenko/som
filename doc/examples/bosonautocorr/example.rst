.. _example_bosonautocorr:

Example: Autocorrelator of a Hermitian operator
===============================================

Observable kind: `BosonAutoCorr`.

Correlator of a Hermitian operator with itself.

For the Hermitian operators the auxiliary spectral function
:math:`A(\epsilon) = \Im\chi(\epsilon)/\epsilon` is non-negative and symmetric.

This mode is faster and more robust than `BosonCorr`.

Run analytical continuation
---------------------------

.. literalinclude:: example.py

Download input file :download:`example.h5`.

Plot input and reconstructed correlators at Matsubara frequencies
-----------------------------------------------------------------

.. plot:: examples/bosonautocorr/plot_chi_iw.py
    :include-source:
    :scale: 100

Plot the correlator on the real frequency axis
----------------------------------------------

.. plot:: examples/bosonautocorr/plot_chi_w.py
    :include-source:
    :scale: 100
