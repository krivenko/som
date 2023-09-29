.. _script_porting:

SOM 1.x script porting guide
============================

.. currentmodule:: som

Porting computational scripts from SOM 1.x to SOM 2.0 means switching from
Python 2.7 to Python 3 and from TRIQS 1.4 to TRIQS 3.1. This guide covers
only the SOM-specific portion of changes that need to be made. Please, refer to
the official Python
`porting guide <https://docs.python.org/3/howto/pyporting.html>`_ to learn about
the general language changes introduced in Python 3, such as ``print()``
becoming a function instead of a statement, and the new semantics of the integer
division. There is also a page about `porting applications
<https://triqs.github.io/triqs/latest/porting_to_triqs3.html>`_ to TRIQS 3.0
provided by TRIQS' developers.

Python modules
~~~~~~~~~~~~~~

Following a convention change for naming of TRIQS applications, the main Python
module of SOM has been renamed.

::

    ### SOM 1.x ###
    from triqs.applications.analytical_continuation.som import Som

::

    ### SOM 2.0 ###
    from som import Som

Functions implementing the :ref:`statistical analysis of ensembles of spectral
functions <spectral_stats>` are collected in a new module,

::

    ### SOM 2.0 ###
    from som.spectral_stats import (spectral_integral,
                                    spectral_avg,
                                    spectral_disp,
                                    spectral_corr)

Construction of the :class:`Som` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- It is now possible to provide full :ref:`covariance matrices <cov_matrix>`
  as an alternative to estimated error bars upon construction of the
  :class:`Som` object.

  ::

      ### SOM 1.x ###
      cont = Som(g, error_bars, kind=kind, norms=norms)

  ::

      ### SOM 2.0 ###
      cont_eb = Som(g, error_bars, kind=kind, norms=norms)
      cont_cm = Som(g,
                    cov_matrices,
                    kind=kind,
                    norms=norms,
                    filtering_levels=1e-5)

  The optional argument ``filtering_levels`` :ref:`improves stability of the
  algorithm <cov_matrix_filtered>` when the covariance matrices are used.

- When continuing fermionic Green's functions, there is an option to enforce
  the :ref:`particle-hole symmetry of the spectrum <fermiongfsymm>` by passing
  ``kind="FermionGfSymm"`` instead of ``kind="FermionGf"``.

- Definition of the integral kernels for ``kind="BosonAutoCorr"`` has been
  changed: The spectral function :math:`A(\epsilon)` is now defined on the
  whole energy axis instead of :math:`\epsilon\in[0;\infty[`, while the kernels
  gained an extra coefficient :math:`1/2`. It has been observed that the new
  definition makes the algorithm better reproduce results of the ``BosonCorr``
  kernels for the same input data. This change has an implication on scripting,
  both ``BosonCorr`` and ``BosonAutoCorr`` expect the same normalization
  constants in ``norms`` from now on (with SOM 1.x one had to divide the
  constants by 2 for ``BosonAutoCorr``).

- If normalization constants for the ``BosonCorr`` or ``BosonAutoCorr``
  spectra are not known a priori, they can be estimated from the input data
  by calling a new utility function :func:`estimate_boson_corr_spectrum_norms`.

  ::

      ### SOM 2.0 ###
      from som import estimate_boson_corr_spectrum_norms

      # Given a correlator of boson-like operators $\chi$ defined on any
      # supported mesh, return a list of spectrum normalization constants
      # $\mathcal{N} = \pi \chi(i\Omega = 0)$.
      norms = estimate_boson_corr_spectrum_norms(chi)

      cont = Som(chi, error_bars, kind="BosonCorr", norms=norms)

Deprecation of :func:`Som.run()`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to accommodate for new features, method :func:`Som.run()` has been
declared deprecated, and its functionality has been split between a few new
methods.

  ::

      ### SOM 1.x ###
      cont.run(**params)

  ::

      ### SOM 2.0 ###

      # Accumulate particular solutions.
      cont.accumulate(**acc_params)

      # Compute the final solution using the procedure from SOM 1.x.
      cont.compute_final_solution(**fs_params)

      # or

      # Compute the final solution using the Consistent Constraints procedure
      # new to SOM 2.0.
      cont.compute_final_solution_cc(**fscc_params)

Passing ``cc_update=True`` to :func:`Som.accumulate` will enable the
:ref:`Consistent constraints update <cc_update>` that may speed up search for
better particular solutions. A bunch of :func:`Som.accumulate`'s parameters
named ``cc_update_*`` give a means to fine-tune behavior of the CC updates.

Calling :func:`Som.accumulate()` multiple times will incrementally extend the
pool of accumulated particular solutions. :func:`Som.clear()` will remove
all accumulated solutions.

In SOM 1.x, :func:`Som.run()` was selecting good particular solutions based on
a criterion established by parameter ``adjust_l_good_d``. With SOM 2.0,
selection of good particular solutions is performed as part of algorithms
implemented in :func:`Som.compute_final_solution` and
:func:`Som.compute_final_solution_cc`.
They both accept arguments ``good_chi_rel`` and ``good_chi_abs``, and
select good solution based on values of the
:ref:`"goodness of fit" <goodness_of_fit>` :math:`\chi^2`-functional associated
with those solutions. A good solution :math:`A_j` must simultaneously satisfy
:math:`\chi[A_j] \leq` ``good_chi_abs`` and
:math:`\chi[A_j] \leq \min_{j'}(\chi[A_{j'}])\times` ``good_chi_rel``.

:func:`Som.compute_final_solution_cc` constructs
the :ref:`final solution <final_solution>` using a sophisticated iterative
optimization procedure with many adjustable parameters. It can result in a
smoother spectral function, which can optionally be biased towards a
user-provided default model.

Automatic adjustment of the number of global updates per solution (:math:`F`),
which used to be one of :func:`Som.run`'s features, is now available as method
:func:`Som.adjust_f`.

::

    ### SOM 2.0 ###

    # Adjust the number of global updates.
    f = cont.adjust_f(energy_window=(-5, 5))
    # Accumulate particular solutions.
    cont.accumulate(energy_window=(-5, 5), f=f, **acc_params)

Post-processing of spectral functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Due to changes in the TRIQS Green's function library, it is no longer possible
to use the ``g << cont`` syntax. Furthermore, information about the high
frequency expansion (tail) coefficients has been separated from Green's function
container objects. The following snippets show the updated syntax for recovering
the real-frequency versions of observables, reconstructing the imaginary
time/Matsubara frequency/Legendre coefficient data and computing the tail.

::

    ### SOM 1.x ###

    # Recover the real-frequency counterpart of 'g' and its tail.
    g_w = GfReFreq(window=energy_window, n_points=n_w, indices=g.indices)
    g_w << cont

    # Reconstruct the input quantity from the computed spectral function.
    g_rec = g.copy()
    g_rec << cont

::

  ### SOM 2.0 ###

  from som import fill_refreq, compute_tail, reconstruct

  # Recover the real-frequency counterpart of 'g'.
  g_w = GfReFreq(window=energy_window, n_points=n_w, indices=g.indices)
  fill_refreq(g_w, cont)

  # Compute the tail of 'g_w'.
  tail = compute_tail(tail_max_order, cont)

  # Reconstruct an observable from the computed spectral function.
  g_rec = g.copy()
  reconstruct(g_rec, cont)

By default, :func:`fill_refreq` uses :ref:`binning <binning>`, which can be
disabled by passing ``with_binning=False``.

Direct access to spectral functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A few new attributes added to :class:`Som` give access to the accumulated
particular solutions, the final solution and their respective values of the
:math:`\chi^2`-functional.

::

    ### SOM 2.0 ###

    # Extract a list of pairs (accumulated particular solution, its \chi^2) for
    # the 0-th diagonal component of the observable.
    # This list is local to the calling MPI rank.
    part_sols_with_chi2 = cont.particular_solutions(0)

    # Minimum of \chi^2 over all accumulated particular solutions
    # on all MPI ranks.
    chi2_min = cont.objf_min

    # List of final solutions, one element per diagonal component of
    # the observable.
    final_sols =  cont.solutions

    # List of \chi^2 values for the final solutions, one element per diagonal
    # component of the observable.
    chi2_final = cont.objf_list

All solutions extracted this way are instances of a new class
:class:`Configuration`, which is a collection of :class:`Rectangle`'s.
Configurations (spectral functions) can be iterated over, evaluated at a given
value of energy, stored to/loaded from an
:ref:`HDF5 archive <triqslibs:hdf5_tutorial>`.
