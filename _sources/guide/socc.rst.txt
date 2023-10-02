.. _socc:

Stochastic Optimization with Consistent Constraints
===================================================

.. currentmodule:: som

Application of the method of consistent constraints in SOM is twofold:

- A new Markov chain update that is enabled by calling :func:`Som.accumulate`
  with ``cc_update=True``.
- A new method :func:`Som.compute_final_solution_cc` to build the
  :ref:`final solution <final_solution>` out of accumulated
  :ref:`particular solutions <particular_solutions>`.

Both applications are rather involved and support many user-adjustable
parameters. Names of their respective keyword arguments are typeset as
``monospaced text`` in the following.

.. _cc_update:

Consistent-constraints updates
------------------------------

Consistent-constraints (CC) updates introduced in SOM 2.0 are a new type of
updates in the Markov chain used to accumulate particular solutions
(Section II.A of [GMPPS2017]_). Unlike elementary updates, they radically change
an MC configuration and help reveal local minima of the objective function more
quickly.

The CC updates are proposed during the first stage of a global update at regular
intervals. Two consecutive CC updates are separated by
``cc_update_cycle_length`` elementary updates. A CC update is skipped if it
could potentially result in a configuration with too many rectangles,
:math:`2K_0+1>` ``max_rects``, where :math:`K_0` is the number of rectangles in
the current configuration. CC updates are accepted and rejected based on the
same Metropolis criterion as the elementary updates. A proposed CC update
proceeds as follows.

#. Take the current configuration, construct a configuration made out of
   non-overlapping rectangles :math:`\{(c_k, w_k, h_k)\}` (Eqs. (16)-(17) of
   [GMPPS2017]_) and discard rectangles with width below ``min_rect_width``.
   If size :math:`K` of the resulting non-overlapping configuration is less than
   3 (too small to evaluate the 2nd derivative), reject the update.
#. Collect heights :math:`h_k` of rectangles in the non-overlapping
   configuration in a :math:`K`-dimensional vector :math:`\mathbf{h}`.
#. Starting from the collected heights, use the CC protocol (see below) to
   compute optimal heights :math:`\mathbf{h}^\mathrm{(opt)}`.
#. If any of :math:`h^\mathrm{(opt)}_k` is significantly negative, i.e.
   :math:`w_k h^\mathrm{(opt)}_k < -` ``min_rect_weight``, reject the update.
#. Replace heights :math:`h_k` in the non-overlapping configuration with the
   optimal heights :math:`h^\mathrm{(opt)}_k`.
#. Remove small rectangles
   (:math:`w_k h^\mathrm{(opt)}_k <` ``min_rect_weight``) from the configuration
   and redistribute their weight by changing heights of their neighboring
   rectangles. This step is repeated until none of the rectangles has weight
   below ``min_rect_weight``.
#. Evaluate the :math:`\chi^2`-functional for the new optimized configuration
   and use :math:`\chi^2_\mathrm{(opt)}` to compute the Metropolis ratio.

The CC protocol consists in iterative minimization of a quadratic function of
:math:`h_k` with self-consistent determination of regularization parameters
:math:`Q^{(0)}_k, Q^{(1)}_k, Q^{(2)}_k` this function depends on. The function
in question is a sum of the :math:`\chi^2`-functional and a few regularization
terms,

.. math::

    O[\mathbf{h}; Q^{(0)},Q^{(1)},Q^{(2)}] = \chi^2[\mathbf{h}]
     + O_0[\mathbf{h};Q^{(0)}]
     + O_1[\mathbf{h};Q^{(1)}]
     + O_2[\mathbf{h};Q^{(2)}].

The regularization term :math:`O_0` penalizes large amplitudes of a solution,

.. math::

    O_0[\mathbf{h};Q^{(0)}] = \sum_{k=1}^K (Q^{(0)}_k)^2 h_k^2 =
        \mathbf{h}^T \hat O_{(0)} \mathbf{h}, \quad
        \hat O_{(0)} = \mathrm{diag}((Q^{(0)}_k)^2).

:math:`O_1` and :math:`O_2` penalize large values of the 1st and 2nd
derivatives of a solution respectively,

.. math::

    O_1[\mathbf{h};Q^{(1)}] &= \sum_{k=1}^{K-1}
        (Q^{(1)}_k)^2 |A'(\epsilon_k)|^2 =
        \mathbf{h}^T \hat O_{(1)} \mathbf{h},

    O_2[\mathbf{h};Q^{(2)}] &= \sum_{k=2}^{K-1}
        (Q^{(2)}_{k-1})^2 |A''(\epsilon_k)|^2 =
        \mathbf{h}^T \hat O_{(2)} \mathbf{h}.

Matrices :math:`\hat O_{(1)}` and :math:`\hat O_{(2)}` are derived from
finite-difference approximations of :math:`A'` and :math:`A''` at points
:math:`\epsilon_k = c_k`.

CC protocol iterations are organized as follows.

#. Initialize regularization parameters according to

   - :math:`Q^{(0)}_k = 0,\ k=\overline{1,K}`,
   - :math:`Q^{(1)}_k = (` ``cc_update_der_penalty_init``
     :math:`)W^2/ \mathcal{N},\ k = \overline{1,K-1}`,
   - :math:`Q^{(2)}_k = (` ``cc_update_der_penalty_init``
     :math:`)W^3/ \mathcal{N},\ k = \overline{1,K-2}`,

   where :math:`W` is the energy window width
   (``energy_window[1] - energy_window[0]``) and :math:`\mathcal{N}` is the
   requested solution norm.

#. Use a numerical linear algebra algorithm to minimize the quadratic function
   :math:`O[\mathbf{h}; Q^{(0)},Q^{(1)},Q^{(2)}]` w.r.t. :math:`\mathbf{h}`
   subject to equality constraint :math:`\sum_{k=1}^K w_k h_k = \mathcal{N}`.

#. Compute discrepancy of weights of rectangles with heights
   :math:`\mathbf{h}` and :math:`\mathbf{h}^\mathrm{new}` as
   :math:`\frac{1}{\mathcal{N}}\sum_{k=1}^K |w_k (h_k - h^\mathrm{new}_k)|`.
   If it lies within a fixed tolerance level
   ``cc_update_rect_norm_variation_tol``, terminate iterations.

#. Update amplitude regularization parameters :math:`Q^{(0)}_k`:
   If :math:`h^\mathrm{new}_k` is negative, set
   :math:`Q^{(0)}_k =(` ``cc_update_height_penalty_max``
   :math:`)W / \mathcal{N}`,
   otherwise divide :math:`Q^{(0)}_k` by ``cc_update_height_penalty_divisor``.

#. Use finite-difference approximations to compute derivatives of :math:`A`
   at :math:`\epsilon_k = c_k`:

   - 1st derivatives :math:`d^{(1)}_k, k = \overline{1,K-1}`.
   - 2nd derivatives :math:`d^{(2)}_k, k = \overline{1,K-2}`.

#. Compute limiting constants for :math:`Q^{(1)}_k` and :math:`Q^{(2)}_k`:

   - :math:`Q^{(1)}_\mathrm{lim} = (`
     ``cc_update_der_penalty_limiter`` :math:`) \min_k Q^{(1)}_k`.
   - :math:`Q^{(2)}_\mathrm{lim} = (`
     ``cc_update_der_penalty_limiter`` :math:`) \min_k Q^{(2)}_k`.

#. Update the derivative regularization parameters based on the magnitudes of
   :math:`d^{(1)}_k` and :math:`d^{(2)}_k`.

   - If :math:`|d^{(1)}_k| > (` ``cc_update_der_penalty_threshold``
     :math:`/ \sqrt{K-1}) / Q^{(1)}_k`, then set
     :math:`Q^{(1)}_k = (` ``cc_update_der_penalty_threshold``
     :math:`/ \sqrt{K-1}) / |d^{(1)}_k|`.
     Otherwise multiply :math:`Q^{(1)}_k` by
     ``cc_update_der_penalty_increase_coeff``.
   - If :math:`|d^{(2)}_k| > (` ``cc_update_der_penalty_threshold``
     :math:`/ \sqrt{K-1}) / Q^{(2)}_k`, then set
     :math:`Q^{(2)}_k = (` ``cc_update_der_penalty_threshold``
     :math:`/ \sqrt{K-1}) / |d^{(2)}_k|`.
     Otherwise multiply :math:`Q^{(2)}_k` by
     ``cc_update_der_penalty_increase_coeff``.

#. Reduce excessively large derivative regularization parameters to avoid
   divergent behaviour:

   - If :math:`Q^{(1)}_k > Q^{(1)}_\mathrm{lim}`, then set
     :math:`Q^{(1)}_k = Q^{(1)}_\mathrm{lim}`.
   - If :math:`Q^{(2)}_k > Q^{(2)}_\mathrm{lim}`, then set
     :math:`Q^{(2)}_k = Q^{(2)}_\mathrm{lim}`.

#. If the maximal allowed number of iterations ``cc_update_max_iter`` is
   reached, terminate. Otherwise set
   :math:`\mathbf{h} = \mathbf{h}^\mathrm{new}` and repeat from step 2.

.. _final_solution_cc:

Construction of final solution using the consistent-constraints protocol
------------------------------------------------------------------------

SOM 2.x features an optimization procedure that finds a vector of
:ref:`final solution <final_solution>` coefficients :math:`\mathbf{c}`
resulting in a smoother spectral function. It also allows to bias the final
solution towards a user-provided default model.

The optimization procedure minimizes a quadratic function of :math:`\mathbf{c}`
depending on a few sets of regularization parameters,

.. math::

  O[\mathbf{c}; Q_k, D_k, B_k, T_k, A_T(\epsilon_k), U] =
        O_Q[\mathbf{c}; Q_k] +
        O_D[\mathbf{c}; D_k] +
        O_B[\mathbf{c}; B_k] +
        O_T[\mathbf{c}; T_k, A_T(\epsilon)] +
        O_U[\mathbf{c}; U].

- :math:`O_Q[\mathbf{c}; Q_k]` stabilizes the final solution by penalizing
  large amplitudes,

  .. math::

    O_Q[\mathbf{c}; Q_k] = \sum_{k=1}^K Q_k^2 A(\epsilon_k)^2 =
    \mathbf{c}^T \hat O_Q \mathbf{c}, \quad
    (\hat O_Q)_{jj'} = \sum_{k=1}^K Q_k^2 A_j(\epsilon_k) A_{j'}(\epsilon_k).

- :math:`O_D[\mathbf{c}; D_k]` penalizes large 1st derivatives of the final
  solution,

  .. math::

      O_D[\mathbf{c}; D_k] = \sum_{k=1}^{K-1} D_k^2 |A'(\epsilon_k)|^2 =
      \mathbf{c}^T \hat O_D \mathbf{c}, \quad
      (\hat O_D)_{jj'} = \sum_{k=1}^{K-1} D_k^2
      A'_j(\epsilon_k) A'_{j'}(\epsilon_k).

- :math:`O_B[\mathbf{c}; B_k]` penalizes large 2nd derivatives of the final
  solution,

  .. math::

      O_B[\mathbf{c}; B_k] = \sum_{k=2}^{K-1} B_{k-1}^2 |A''(\epsilon_k)|^2 =
      \mathbf{c}^T \hat O_B \mathbf{c}, \quad
      (\hat O_B)_{jj'} = \sum_{k=2}^{K-1} B_{k-1}^2
      A''_j(\epsilon_k) A''_{j'}(\epsilon_k).

- :math:`O_T[\mathbf{c}; T_k, A_T(\epsilon)]` penalizes deviations from a
  predefined default model ("target") :math:`A_T(\epsilon)`,

  .. math::

    O_T[\mathbf{c}; T_k, A_T(\epsilon)] &= \sum_{k=1}^K T_k
    [A(\epsilon_k) - A_T(\epsilon_k)]^2 =
    \mathbf{c}^T \hat O_T \mathbf{c} - 2\mathbf{f}_T^T\mathbf{c}+\mathrm{const},

    (\hat O_T)_{jj'} &= \sum_{k=1}^K T_k A_j(\epsilon_k) A_{j'}(\epsilon_k),
    \quad
    (\mathbf{f}_T)_j = \sum_{k=1}^K T_k A_j(\epsilon_k) A_T(\epsilon_k).

- :math:`O_U[\mathbf{c}; U]` penalizes deviations from the equal weight
  superposition :math:`c_j = 1/J`.

  .. math::

      O_U[\mathbf{c}; U] &= U \sum_{j=1}^J (c_j - 1/J)^2 =
      \mathbf{c}^T \hat O_U \mathbf{c}-2\mathbf{f}_U^T\mathbf{c}+\mathrm{const},

      (\hat O_U)_{jj'} &= U \delta_{jj'},\quad (\mathbf{f}_U)_j = U / J.

Values of particular spectral functions :math:`A_j(\epsilon)` and their
derivatives are approximated on a user-supplied uniform or non-uniform energy
grid :math:`\{\epsilon_k\}` (``refreq_mesh``).

The iterative optimization procedure proceeds as follows.

#. Precompute approximate values of :math:`A_j(\epsilon_k)`,
   :math:`A'_j(\epsilon_k)` and :math:`A''_j(\epsilon_k)`.

#. Using user-supplied parameters :math:`A_T(\epsilon_k)=` ``default_model``,
   :math:`T_k =` ``default_model_weights`` and :math:`U=` ``ew_penalty_coeff``
   precompute :math:`\mathbf{f} = \mathbf{f}_T + \mathbf{f}_U`.

#. Initialize regularization parameters according to

   - :math:`\mathcal{D} = (` ``der_penalty_init`` :math:`) / J`,
   - :math:`Q_k = 0`,
   - :math:`D_k = \mathcal{D}`,
   - :math:`B_k = \mathcal{D}`.

#. Use the regularization parameters to compute matrix :math:`\hat O =
   \hat O_Q + \hat O_D + \hat O_B + \hat O_T + \hat O_U`.

#. Use a numerical linear algebra algorithm to minimize
   :math:`\mathbf{c}^T \hat O \mathbf{c} - 2\mathbf{f}^T\mathbf{c}` w.r.t.
   :math:`\mathbf{c}` subject to the normalization sum rule
   :math:`\sum_{J=1}^J c_j = 1`. The minimization result at the :math:`I`-th
   iteration is stored as a vector :math:`\mathbf{c}^{(I)}`.

#. Compute and store relative :math:`L_1`-distance between
   :math:`\mathbf{c}^{(I)}` and :math:`\mathbf{c}^{(I-1)}`,

   .. math::

      d(I, I-1) = \frac{\sum_{j=1}^J |c^{(I)}_j - c^{(I-1)}_j|}
                       {\sum_{j=1}^J |c^{(I)}_j|},

   where the initial vector :math:`\mathbf{c}^{(0)}` is taken to have all
   components equal to :math:`1/J`.

#. In case user has provided a monitor function (``monitor``), it is called
   with 4 arguments,

   - a 1D NumPy array :math:`\mathbf{c}^{(I)}`,
   - a 2-tuple of 1D arrays :math:`(A(\epsilon_k), Q_k)`,
   - a 2-tuple of 1D arrays :math:`(A'(\epsilon_k), D_k)`,
   - a 2-tuple of 1D arrays :math:`(A''(\epsilon_k), B_k)`.

   Returning ``True`` from the function terminates iterations.

#. If :math:`\sum_{j=1}^J |c^{(I)}_j| >` ``max_sum_abs_c`` (contribution of the
   negative coefficients is too big), then terminate iterations.

#. If :math:`d(I, I-1) < \min_n |\sigma_n/G(n)|`, then convergence is reached
   and iterations are terminated. :math:`G(n)` is the right-hand side of the
   analytic continuation problem, and :math:`\sigma_n` are respective
   estimated error bars on the input data.

#. Form a new final solution from coefficients :math:`\mathbf{c}^{(I)}`,

.. math::

      A(\epsilon_k) = \sum_{j=1}^J c^{(I)}_j A_j(\epsilon_k),\quad
      A'(\epsilon_k) = \sum_{j=1}^J c^{(I)}_j A'_j(\epsilon_k),\quad
      A''(\epsilon_k) = \sum_{j=1}^J c^{(I)}_j A''_j(\epsilon_k).

#. Increase regularization parameters :math:`\mathcal{D}`, :math:`D_k` and
   :math:`B_k` by the factor ``der_penalty_coeff``.

#. Update the amplitude regularization parameters :math:`Q_k`: If
   :math:`A(\epsilon_k) < 0`, then set :math:`Q_k =` ``amp_penalty_max``.
   Otherwise divide :math:`Q_k` by ``amp_penalty_divisor``.

#. Update the derivative regularization parameters:

   - If :math:`D_k |A'(\epsilon_k)| > \mathcal{D}`, then set
     :math:`D_k = \mathcal{D} / |A'(\epsilon_k)|`;
   - If :math:`B_k |A''(\epsilon_k)| > \mathcal{D}`, then set
     :math:`B_k = \mathcal{D} / |A''(\epsilon_k)|`.

#. If the maximal allowed number of iterations ``max_iter`` is reached,
   terminate. Otherwise repeat from step 4.

Once the iterations have been terminated, one inspects the accumulated list of
pairs :math:`(\mathbf{c}^{(1)}, d(1,0)), (\mathbf{c}^{(2)}, d(2,1)), \ldots`.
Elements of the list are sorted in the ascending order according to
:math:`d(I,I-1)`, i.e. :math:`\mathbf{c}^{(I)}` from the iteration closest to
convergence comes first in the sorted list. The first vector of coefficients
from the sorted list that obeys

- :math:`\chi^2\left[\sum_{j=1}^J c^{(I)}_j A_j\right] \leq (`
  ``good_chi_abs`` :math:`)^2` and
- :math:`\chi^2\left[\sum_{j=1}^J c^{(I)}_j A_j\right] \leq \min_j \chi^2[A_j]
  (` ``good_chi_rel`` :math:`)^2`

represents the sought final solution.
