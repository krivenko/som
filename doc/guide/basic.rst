.. _basic:

Basic description of the method
===============================

.. _integral_equation:

Integral equation of the analytic continuation problem
------------------------------------------------------

The analytic continuation problem amounts to solving a Fredholm integral
equation of the first kind with an approximately known right hand side
:math:`G(n)`,

.. math::

  \int_{\epsilon_\mathrm{min}}^{\epsilon_\mathrm{max}}
  K(n, \epsilon) A(\epsilon) d\epsilon = G(n).

The integral kernel :math:`K(n, \epsilon)` and the integration interval
:math:`[\epsilon_\mathrm{min}; \epsilon_\mathrm{max}]` depend on what
:ref:`kind of physical observable <observables>` :math:`G(n)` is (Green's
function, susceptibility, conductivity). The formal discrete index
:math:`n = \overline{1, N}` may denote one of the following variables,

- Fermionic of bosonic
  :class:`Matsubara frequencies <triqslibs:triqs.gf.meshes.MeshImFreq>`;
- Points of an :class:`imaginary time grid <triqs.gf.meshes.MeshImTime>`;
- Indices of
  :class:`Legendre orthogonal polynomials <triqs.gf.meshes.MeshLegendre>`.

.. _solution_norm:

The integral equation is solved with respect to a non-negative spectral function
:math:`A(\epsilon)` subject to a normalization condition

.. math::

  \int_{\epsilon_\mathrm{min}}^{\epsilon_\mathrm{max}}
  A(\epsilon) d\epsilon = \mathcal{N}.

Goodness of fit functionals
---------------------------

The right hand side :math:`G(n)` is typically obtained via statistical
(Monte Carlo) sampling and averaging of a random quantity,

.. math::

  G(n) = \frac{1}{S}\sum_{s=1}^S G^{(s)}(n).

In addition to this average, SOM requires knowledge of either estimated error
bars

.. math::

  \sigma(n) = \sqrt{\frac{1}{S(S-1)} \sum_{s=1}^S |G^{(s)}(n) - G(n)|^2}

or of a full covariance matrix

.. math::

  \Sigma_{nn'} = \frac{1}{S(S-1)} \sum_{s=1}^S (G^{(s)}(n) - G(n))^*
                                               (G^{(s)}(n') - G(n')).

.. _error_bars:

When only the error bars are known, SOM tries to minimize the "goodness of fit"
functional

.. math::
  \chi^2[A(\epsilon)] = \frac{1}{N} \sum_{n=1}^N
  \left|\frac{\Delta(n)}{\sigma(n)}\right|^2,

where :math:`\Delta(n)` are discrepancies

.. math::
    \Delta(n) = \int_{\epsilon_\mathrm{min}}^{\epsilon_\mathrm{max}}
                K(n, \epsilon) A(\epsilon) d\epsilon - G(n).

.. _cov_matrix:

With the covariance matrix available, a more general functional can be used
instead,

.. math::

    \chi^2[A(\epsilon)] = \frac{1}{N}
    \mathbf{\Delta}^\dagger \hat\Sigma^{-1} \mathbf{\Delta}, \quad
    \mathbf{\Delta} = \{\Delta(n)\}.

.. _cov_matrix_filtered:

Due to statistical noise, some of :math:`\hat\Sigma`'s eigenvalues may turn
to be negative, which would make :math:`\chi^2` lack a global minimum.
Furthermore, very small positive eigenvalues are still problematic as they can
lead to numerical instability of the minimization procedure.
SOM tackles these problems by discarding all negative eigenvalues and shifting
the positive ones by a constant
:math:`l^2` (:math:`l` is referred to as "filtering level"). This
results in a modified :math:`\chi^2`-functional

.. math::

    \chi^2[A(\epsilon)] = \frac{1}{\tilde N}
    \mathbf{\Delta}^\dagger \hat{\mathfrak{S}}^{-1} \mathbf{\Delta},

where :math:`\tilde N` is the number of the retained positive eigenvalues
:math:`\sigma^2(n)` and :math:`\hat{\mathfrak{S}}^{-1}` is a regularized version
of :math:`\hat\Sigma^{-1}`,

.. math::

    \hat{\mathfrak{S}}^{-1} = \hat U \mathrm{diag}
    \left(\frac{1}{\sigma^2(n) + l^2}\right)
    \hat U^\dagger.

Here, :math:`\hat U` is a :math:`N{\times}\tilde N` matrix, whose columns are
eigenvectors of :math:`\hat{\Sigma}` corresponding to the retained eigenvalues.

Minimization of the :math:`\chi^2`-functionals
----------------------------------------------

.. _particular_solutions:

The :math:`\chi^2`-minimization proceeds in two steps.

At first, a stochastic minimization algorithm is used to accumulate :math:`J`
*particular solutions* corresponding to local minima of :math:`\chi^2`.
Each solution is parameterized as a superposition of rectangles with positive
heights,

.. math::

    A(\epsilon) = \sum_{k=1}^K R_{\{c_k, w_k, h_k\}}(\epsilon),

    R_{\{c, w, h\}}(\epsilon) \equiv h \theta(w/2-|\epsilon-c|).

Accumulation of each solution is performed using a Markov chain consisting of
:math:`F` global updates. In turn, each global update comprises :math:`T`
elementary updates that continuously change center positions :math:`c_k`,
widths :math:`w_k`, heights :math:`h_k` as well the total number of rectangles
:math:`K`. Further acceleration of the optimization procedure can be achieved
by enabling the large scale :ref:`Consistent Constraints updates <cc_update>`.

.. _final_solution:

At the second step, :math:`\tilde J` out of the :math:`J` particular solutions
are classified as good based on their :math:`\chi^2[A_j(\epsilon)]`: A good
solution must satisfy :math:`\chi[A_j] \leq \chi_c^\mathrm{abs}` and
:math:`\chi[A_j] \leq \min_j(\chi[A_j]) \chi_c^\mathrm{rel}`.

The good solutions are then combined to form a *final solution*,

.. math::

  A(\epsilon) = \sum_{j=1}^{\tilde J} c_j A_j(\epsilon), \quad
  \sum_{j=1}^{\tilde J} c_j = 1

In the traditional formulation of SOM [MPSS2000]_, the particular solutions are
summed up with equal weights, :math:`c_j = 1/\tilde J`. It is also possible to
employ a more sophisticated and customizable :ref:`iterative procedure
<final_solution_cc>` proposed in [GMPPS2017]_. It yields a set of coefficients
:math:`c_j` that make for smoother spectra, and can optionally favor solutions
close to a user-defined default model.

Instead of computing the final solution, one can use the
:ref:`statistical analysis technique <spectral_stats>` from Sections I-II of
[GMPPS2017]_ to estimate averages, dispersions and 2-point correlations of
the particular solution ensemble.

Post-processing of spectral functions
-------------------------------------

There are a few operations one can perform with a computed final spectral
function :math:`A(\epsilon)`.

- :func:`Extract <som.Som.solutions>` and :class:`inspect individual rectangles
  <som.Configuration>` of the final solutions in their sum-of-rectangles form.

- :ref:`Recover the real-frequency version <recovery>` of the studied
  observable by projecting it onto a :class:`real frequency mesh
  <triqslibs:triqs.gf.meshes.MeshReFreq>`.
  For instance, the retarded Green's function :math:`G^\mathcal{ret}(\epsilon)`
  can be recovered from a fermionic Matsubara Green's function :math:`G(\tau)`.

- :ref:`Compute the high frequency expansion coefficients <compute_tail>`
  of the real-frequency observable.

- Back-substitute :math:`A(\epsilon)` into the original integral equation and
  reconstruct the rights hand side (not necessarily on the same :math:`n`-mesh).
  This feature is useful for a quick assessment of quality of the solution.
