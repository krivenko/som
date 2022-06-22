.. _example_fermiongf:

Fermionic Green's function or self-energy
=========================================

.. note::

    One can continue a self-energy function as long as it does not contain
    a static Hartree-Fock contribution (i.e. decays to 0 at
    :math:`\omega\to\infty`). In this case norms must be precomputed separately
    as first spectral moments of the self-energy.

    For derivation of spectral moments see, for instance,

    .. code-block:: text

        "Interpolating self-energy of the infinite-dimensional Hubbard model:
        Modifying the iterative perturbation theory"
        M. Potthoff, T. Wegner, and W. Nolting, Phys. Rev. B 55, 16132 (1997)

Perform analytic continuation using the :ref:`FermionGf <fermiongf>` kernel
---------------------------------------------------------------------------

.. literalinclude:: example.py

Download input file :download:`example.h5`.

Plot input and reconstructed imaginary-time GF's
------------------------------------------------

.. plot:: examples/fermiongf/plot_g_tau.py
    :include-source:
    :scale: 100

Plot spectral functions and tail coefficients
---------------------------------------------

.. plot:: examples/fermiongf/plot_g_w.py
    :include-source:
    :scale: 100

Plot :math:`\chi`-histograms
----------------------------

.. plot:: examples/fermiongf/plot_hist.py
    :include-source:
    :scale: 100
