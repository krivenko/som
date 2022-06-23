.. _example_zerotemp:

Dynamical response function at zero temperature
===============================================

Formally speaking, the imaginary time segment :math:`\tau\in[0;\beta)`
turns into an infinite interval as :math:`\beta\to\infty`. Similarly,
spacing between Matsubara frequencies goes to 0 in that limit, and
the difference between fermionic and bosonic Matsubaras disappears.

One can still define a correlation function on a finite time mesh
:math:`\tau_i\in[0;\tau_{max}]`, and assume the function is zero
for :math:`\tau>\tau_{max}`. In the frequency representation this
corresponds to fictitious Matsubara spacing :math:`2\pi/\tau_{max}`.

The spectral function is defined only on the positive half-axis of energy,
since :math:`(1\pm e^{-\beta\epsilon})^{-1}` vanishes for negative
:math:`\epsilon` in the zero temperature limit.

Perform analytic continuation using the :ref:`ZeroTemp <zerotemp>` kernel
-------------------------------------------------------------------------

.. literalinclude:: example.py

Download input file :download:`example.h5`.

Plot input and reconstructed imaginary-time correlators
-------------------------------------------------------

.. plot:: examples/zerotemp/plot_g_tau.py
    :include-source:
    :scale: 100

Plot spectral function and tail coefficients
--------------------------------------------

.. plot:: examples/zerotemp/plot_g_w.py
    :include-source:
    :scale: 100

Plot :math:`\chi`-histogram
---------------------------

.. plot:: examples/zerotemp/plot_hist.py
    :include-source:
    :scale: 100
