.. _example_meshes:

Input data defined on various meshes
====================================

In this example we analytically continue the same dynamical susceptibility
:math:`\chi`, whose values are given on three different meshes:

- Bosonic :class:`Matsubara frequencies <triqs.mesh.MeshImFreq>`;
- Points of an :class:`imaginary time grid <triqs.mesh.MeshImTime>`;
- Indices of
  :class:`Legendre orthogonal polynomials <triqs.mesh.MeshLegendre>`.

.. literalinclude:: example.py

Download input file :download:`input.h5`.

.. rubric:: Plot input and reconstructed susceptibility

.. plot:: examples/meshes/plot_chi_rec.py
    :include-source:
    :scale: 100

.. rubric:: Plot spectral functions

.. plot:: examples/meshes/plot_chi_w.py
    :include-source:
    :scale: 100
