.. _welcome:

Stochastic Optimization Method for Analytic Continuation
========================================================

SOM is a :ref:`TRIQS-based <triqslibs:welcome>` implementation of the Stochastic
Optimization Method for Analytic Continuation [MPSS2000]_, which solves a family
of Fredholm integral equations of the first kind. Numerical solution of such
equations is required to reconstruct a real-frequency spectral representation of
physical observables (Green's functions, dynamical susceptibilities) from noisy
results of Quantum Monte Carlo simulations.

Said integral equations have the following general form,

.. math::
    \int_{\epsilon_\mathrm{min}}^{\epsilon_\mathrm{max}}
    K(n, \epsilon) A(\epsilon) d\epsilon = G(n),

where :math:`G(n)` is an input quantity, known with some uncertainty and
:math:`A(\epsilon)` is the spectral function to be found. Meaning of the index
:math:`n`, limits of the integral, and the form of the kernel
:math:`K(n, \epsilon)` depend on the problem at hand.
:math:`A(\epsilon)` is assumed to be non-negative and normalized to a known
constant :math:`\mathcal{N}`. A full list of supported observables with
respective integral kernels is given :ref:`here <kernels>`.
All equations of the considered family share a common property of being
*ill-posed problems*: Their solutions are not unique  and tiny variations of the
right hand part may result in huge changes of the solution.

The Stochastic Optimization Method uses Markov chains combined with other
optimization techniques to minimize one of the following objective functions
("goodness of fit" functionals),

.. math::
    \chi^2_\sigma[A(\epsilon)] = \frac{1}{N} \sum_{n=1}^N
    \left|\frac{\Delta(n)}{\sigma(n)}\right|^2, \quad
    \chi^2_\Sigma[A(\epsilon)] = \frac{1}{N}
    \mathbf{\Delta}^\dagger \hat\Sigma^{-1} \mathbf{\Delta},

where :math:`\mathbf{\Delta} = \{\Delta(n)\}` are discrepancies

.. math::
    \Delta(n) = \int_{\epsilon_\mathrm{min}}^{\epsilon_\mathrm{max}}
                K(n, \epsilon) A(\epsilon) d\epsilon - G(n).

These two :math:`\chi^2`-functionals correspond to the two cases when either the
estimated error bars :math:`\sigma(n)` or the :math:`N\times N` covariance
matrix :math:`\hat\Sigma` of the input data :math:`G(n)` are known.

An important feature of the algorithm is its "continuous-energy" representation
of solutions. Function :math:`A(\epsilon)` is parameterized as a sum of
rectangles with certain positions, widths and heights.

SOM release 2.0.0 brings a number of new features, most notably a complete
implementation of the Stochastic Optimization with Consistent Constraints (SOCC)
protocol proposed in [GMPPS2017]_.

.. note::

    This code uses MPI parallelism and is considerably more CPU-intensive than
    algorithms based on the Bayesian statistical inference (MaxEnt) or Padé
    approximants. Its primary design goal is the ability to resolve fine
    spectral features, rather than giving a quick insight into general shape of
    the spectrum.

.. toctree::
    :name: mastertoc
    :maxdepth: 3
    :hidden:

    install
    tutorial
    guide/index
    examples
    reference
    script_porting
    issues
    about
    ChangeLog.md
    genindex
    search

.. [MPSS2000]
   | "Diagrammatic quantum Monte Carlo study of the Fröhlich polaron",
   | A. S. Mishchenko, N. V. Prokof'ev, A. Sakamoto, and B. V. Svistunov,
   | Phys. Rev. B **62**, 6317 (2000)
   | https://doi.org/10.1103/PhysRevB.62.6317

.. [GMPPS2017]
   | "Numerical analytic continuation: Answers to well-posed questions",
   | O. Goulko, A. S. Mishchenko, L. Pollet, N. Prokof'ev, and B. Svistunov,
   | Phys. Rev. B **95**, 014102 (2017)
   | https://doi.org/10.1103/PhysRevB.95.014102

