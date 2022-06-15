.. _socc:

Stochastic Optimization with Consistent Constraints
===================================================

.. currentmodule:: som

Application of the method of consistent constraints in SOM is twofold:

- A new Markov chain update that is enabled by calling :func:`Som.accumulate`
  with ``cc_update=True``.
- A new method :func:`Som.compute_final_solution_cc` to build a
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

#. Take the current configuration and construct a configuration made out of
   non-overlapping rectangles :math:`\{(c_k, w_k, h_k)\}` (Eqs. (16)-(17) of
   [GMPPS2017]_) and discard rectangles with width below ``min_rect_width``.
   If size :math:`K` of the resulting non-overlapping configuration is less than
   3, reject the update.
#. Collect heights :math:`h_k` of rectangles in the non-overlapping
   configuration in a :math:`K` dimensional vector :math:`\mathbf{h}`.
#. Starting from the collected heights, use the CC protocol (see below) to
   compute optimal heights :math:`\mathbf{h}^\mathrm{(opt)}`.
#. If any of :math:`h^\mathrm{(opt)}_k` is significantly negative, i.e.
   :math:`w_k h^\mathrm{(opt)}_k < -` ``min_rect_weight``, reject the update.
#. Replace heights :math:`h_k` in the non-overlapping configuration with the
   optimal heights :math:`h^\mathrm{(opt)}_k`.
#. Remove small rectangles
   (:math:`w_k h^\mathrm{(opt)}_k <` ``min_rect_weight``) from the configuration
   and redistribute their weight by changing heights of their neighboring
   rectangles. This step is repeated until all rectangles satisfy
   (:math:`w_k h_k \geq` ``min_rect_weight``).
#. Evaluate the :math:`\chi^2`-functional for the new optimized configuration,
   and use :math:`Z(\chi^2 /\chi^2_\mathrm{(opt)})` as the Metropolis ratio.

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

:math:`O_1` and :math:`O_2` penalize large values of the 1st and the 2nd
derivative of a solution respectively,

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
