.. _example_spectral_stats:

Statistical analysis of ensembles of spectral functions
=======================================================

This example is an illustration of the :ref:`spectral statistical analysis API
<spectral_stats>`. In the script below, we use the same input data for a
two-component fermionic Green's function as in
:ref:`another example <example_fermiongf>`, but never construct a final
solution out of accumulated particular solutions. Instead, we estimate it by
computing spectral averages over the ensemble of the particular solutions, along
with error bars (dispersions) and two-point correlations.
Use of all three :ref:`resolution functions <resolution_functions>` is shown.

.. rubric:: Accumulate particular solutions and compute
            statistical characteristics

.. literalinclude:: example.py

Download input file :download:`input.h5`.

.. rubric:: Plot statistical characteristics

.. plot:: examples/spectral_stats/plot.py
    :include-source:
    :scale: 100
