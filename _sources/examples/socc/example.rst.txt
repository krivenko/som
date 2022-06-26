.. _example_socc:

Stochastic Optimization with Consistent Constraints
===================================================

Both aspects of SOCC functionality are covered by this example: It shows how to
enable the :ref:`consistent-constraints updates <cc_update>` and how to invoke
the :ref:`SOCC self-consistent optimization algorithm <final_solution_cc>`
to construct the final solution.

When calculating the final solution, we use a predefined default model
("target") and a monitor function that stores SOCC regularization parameters
as they change from one iteration of the optimization process to the next.

.. rubric:: Accumulate particular solutions and apply both SOM and SOCC
            procedures to construct the final solution

.. literalinclude:: example.py

.. rubric:: Plot spectral functions, reconstructed Green's functions and
            evolution of SOCC regularization parameters

.. plot:: examples/socc/plot.py
    :include-source:
    :scale: 100
