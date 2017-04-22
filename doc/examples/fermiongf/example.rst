.. _example_fermiongf:

Example: Fermionic Green's function or self-energy
==================================================

Observable kind: `FermionGf`.

.. note::

    With `FermionGf` one can also continue a self-energy function as long as it
    does not contain a static Hartree-Fock contribution (decays to 0 at :math:`\omega\to\infty`).
    In this case norms must be precomputed separately as first spectral moments of
    the self-energy.

    For derivation of spectral moments see, for instance,

    ::

        "Interpolating self-energy of the infinite-dimensional Hubbard model:
        Modifying the iterative perturbation theory"
        M. Potthoff, T. Wegner, and W. Nolting, Phys. Rev. B 55, 16132 (1997)


Run analytical continuation
---------------------------

.. literalinclude:: example.py

Download input file :download:`example.h5`.

Plot input and reconstructed imaginary-time GF's
------------------------------------------------

.. plot:: examples/fermiongf/plot_g_tau.py
    :include-source:
    :scale: 100

Plot the spectral function
--------------------------

.. plot:: examples/fermiongf/plot_g_w.py
    :include-source:
    :scale: 100
