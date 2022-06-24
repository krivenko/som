.. _example_cov_matrix:

Full covariance matrix of input data
====================================

This example demonstrates three ways to provide information about uncertainty
of input data (here, of an imaginary time fermionic Green's function
:math:`G(\tau)`),

- As :ref:`estimated error bars <error_bars>` :math:`\sigma(\tau)`;
- As a :ref:`full covariance matrix <cov_matrix>` :math:`\hat\Sigma_{\tau\tau'}`
  ;
- As a full covariance matrix :math:`\hat\Sigma_{\tau\tau'}` with all
  eigenvalues shifted up by a constant :math:`l^2`, where :math:`l` is so called
  :ref:`filtering level <cov_matrix_filtered>`.

Perform analytic continuation
-----------------------------

.. literalinclude:: example.py

Download input file :download:`input.h5`.

Plot input and reconstructed imaginary-time GF's
------------------------------------------------

.. plot:: examples/cov_matrix/plot_g_tau.py
    :include-source:
    :scale: 100

Plot spectral functions
-----------------------

.. plot:: examples/cov_matrix/plot_g_w.py
    :include-source:
    :scale: 100
