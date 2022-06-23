.. _example_fermiongfsymm:

Fermionic Green's function or self-energy with enforced particle-hole symmetry
==============================================================================

This example is very similar to :ref:`the previous one <example_fermiongf>`.
This time, however, we choose a special observable kind (``FermionGfSymm``)
that enforces the particle-hole symmetry of the resulting spectral function.

Perform analytic continuation using the :ref:`FermionGfSymm <fermiongfsymm>` kernel
-----------------------------------------------------------------------------------

.. literalinclude:: example.py

Download input file :download:`example.h5`.

Plot input and reconstructed imaginary-time GF's
------------------------------------------------

.. plot:: examples/fermiongfsymm/plot_g_tau.py
    :include-source:
    :scale: 100

Plot spectral function and tail coefficients
--------------------------------------------

.. plot:: examples/fermiongfsymm/plot_g_w.py
    :include-source:
    :scale: 100

Plot :math:`\chi`-histogram
---------------------------

.. plot:: examples/fermiongfsymm/plot_hist.py
    :include-source:
    :scale: 100
