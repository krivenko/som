.. _welcome:

Stochastic Optimization Method for Analytic Continuation
========================================================

The :ref:`TRIQS-based <triqslibs:welcome>` Stochastic Optimization Method for Analytic Continuation allows to solve
a family of Fredholm integral equations of the first kind. Numerical solution of such equations is
required to reconstruct the spectral representation of physical observables (Green's functions,
dynamical susceptibilities) from the noisy results of Quantum Monte Carlo simulations.

The integral equations have the following general form,

.. math::
    \int_{\epsilon_\mathrm{min}}^{\epsilon_\mathrm{max}}
    K(m,\epsilon) A(\epsilon) d\epsilon = G(m),

where :math:`G(m)` is an input quantity, known with some uncertainty, and
:math:`A(\epsilon)` is a spectral function to be found. Meaning of the index :math:`m`,
limits of the integral, and the form of the kernel :math:`K(m,\epsilon)` depend on the problem at hand.
:math:`A(\epsilon)` is assumed to be non-negative and normalized to a known constant.
A full list of supported observables with respective integral kernels is given :ref:`here <kernels>`.

The whole family of equations shares a common property of being *ill-posed problems*;
tiny variations of the right hand part may result in huge changes of the solution at high energies.

The Stochastic Optimization Method uses Markov chain sampling to minimize an objective function

.. math::
    \mathcal{D} = \sum_m \frac{1}{S(m)}\left|
        \int_{\epsilon_\mathrm{min}}^{\epsilon_\mathrm{max}}
        K(m,\epsilon) A(\epsilon) d\epsilon - G(m)
    \right|.

One can use a weight factor :math:`S(m)` to stress importance of some data points :math:`G(m)`.

An important feature of the algorithm is a "continuous-energy" representation of the solutions.
Function :math:`A(\epsilon)` is parametrized as a sum of rectangles with certain positions, widths and
heights.

*Despite a number of optimizations, this code is considerably more CPU-intensive than algorithms based
on the Bayesian statistical inference (MaxEnt) or Padé approximants. Its primary design goal is
the ability to resolve fine spectral features, rather than giving a quick insight into general shape
of the spectrum.*
