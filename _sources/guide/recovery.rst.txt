.. _recovery:

Recovery of observables on a real-frequency mesh
================================================

Recovery of the real-frequency version of a dynamical observable, such as the
:ref:`retarded Green's function <fermiongf>` :math:`G^\mathrm{ret}(\epsilon)`,
from a spectral function :math:`A(\epsilon)` is performed by function
:func:`som.fill_refreq`. It applies a variant of the Hilbert transform to the
spectrum and projects the result onto a regular grid of real energy points

.. math::

  \epsilon_i = \epsilon_\mathrm{min} + i\Delta\epsilon,\
  i=\overline{0, I-1},\
  \Delta\epsilon = \frac{\epsilon_\mathrm{max} - \epsilon_\mathrm{min}}{I - 1}.

A final solution produced by SOM has a general sum-of-rectangles form

.. math::

    A(\epsilon) = \sum_{k=1}^K R_{\{c_k, w_k, h_k\}}(\epsilon),

    R_{\{c, w, h\}}(\epsilon) \equiv h \theta(w/2-|\epsilon-c|).

Before SOM 2.0, the retarded Green's function would be recovered and projected
onto the frequency mesh as

.. math::

    G^\mathrm{ret}(\epsilon_i) = -\int\limits_{-\infty}^\infty d\epsilon'
      \frac{A(\epsilon')}{\epsilon' - \epsilon_i - i0} =
      -\sum_{k=1}^K h_k \int\limits_{c_k-w_k/2}^{c_k+w_k/2}
      \frac{d\epsilon'}{\epsilon' - \epsilon_i - i0}.

.. _binning:

By default, SOM 2.x uses binning while projecting onto the frequency mesh.
It splits the energy window
:math:`[\epsilon_\mathrm{min};\epsilon_\mathrm{max}]` into :math:`I` bins,

.. math::

  \mathfrak{B}_i = \left\{
  \begin{array}{ll}
      [\epsilon_\mathrm{min}; \epsilon_\mathrm{min}+\Delta\epsilon/2], &i=0,\\
      [\epsilon_i-\Delta\epsilon/2; \epsilon_i+\Delta\epsilon/2],
        &i=\overline{1,I-2},\\
      [\epsilon_\mathrm{max}-\Delta\epsilon/2;
       \epsilon_\mathrm{max}+\Delta\epsilon/2], &i=I-1.
  \end{array}
  \right.

and averages the integration results over the bins,

.. math::

    G^\mathrm{ret}(\epsilon_i) \approx
    - \frac{1}{|\mathfrak{B}_i|}
      \int_{\mathfrak{B}_i} d\epsilon
      \int\limits_{-\infty}^\infty d\epsilon'
      \frac{A(\epsilon')}{\epsilon' - \epsilon - i0} =\\=
      -\sum_{k=1}^K h_k
      \frac{1}{|\mathfrak{B}_i|}
      \int_{\mathfrak{B}_i} d\epsilon
      \int\limits_{c_k-w_k/2}^{c_k+w_k/2}
      \frac{d\epsilon'}{\epsilon' - \epsilon - i0}.

Owing to the extra integral, projection with binning makes for smoother
resulting curves for the same number of energy points :math:`\epsilon_i`.
The old behavior can be enabled by passing ``with_binning=False`` to
:func:`som.fill_refreq`.

.. _compute_tail:

:func:`som.compute_tail` is another useful function that extracts information
about the :ref:`high frequency expansion coefficients <triqslibs:gf_tail>`
of :math:`G^\mathrm{ret}(\epsilon)` from :math:`A(\epsilon)`.
