.. _about:

Authors
=======

This analytic continuation program is written by Igor Krivenko (`accompanying paper <https://doi.org/10.1016/j.cpc.2019.01.021>`_,
`arXiv:1808.00603 <https://arxiv.org/abs/1808.00603>`_).

It is based on a method devised by Andrey S. Mishchenko and coauthors.

    `A. S. Mishchenko, N. V. Prokof’ev, A. Sakamoto, and B. V. Svistunov, Phys. Rev. B 62, 6317 <http://dx.doi.org/10.1103/PhysRevB.62.6317>`_

You can find a detailed description of the method in Chapter 14 of the Lecture Notes
from the Autumn School "Correlated Electrons: From Models to Materials".

::

    Eva Pavarini, Erik Koch, Frithjof Anders, and Mark Jarrell (eds.)
    Correlated Electrons: From Models to Materials
    Modeling and Simulation, Vol. 2
    Verlag des Forschungszentrum Jülich, 2012
    ISBN 978-3-89336-796-2

The Lecture Notes are available free of charge online,

    http://www.cond-mat.de/events/correl12/manuscripts/

Acknowledgements
================

I am very grateful to Malte Harland for extensive and thorough testing of the code
at late development stages.

I would also like to acknowledge support from the University of Hamburg, as well as
from the Deutsche Forschungsgemeinschaft via the Project SFB-668
"Magnetismus vom Einzelatom zur Nanostruktur" (A3). Some large-scale tests have been
run on the JURECA HPC machine of the Forschungszentrum Jülich
(project "Continuous Time Quantum Monte Carlo for materials").

TRIQS
=====

It is highly recommended to get familiar with basic usage of the TRIQS library,
before starting to work with this program.

https://triqs.github.io/triqs/1.4/

Here is a simple tutorial showing how to manipulate Green's Function objects of TRIQS
in Python scripts,

https://triqs.github.io/triqs/1.4/tutorials/gfs/gfs_tutorial_python.html

License
=======

The SOM program is published under the `GNU General Public License, version 3
<http://www.gnu.org/licenses/gpl.html>`_.

Note that it *implies* that applications using SOM must also be GPL.

Disclaimer
==========

The program is provided as is, i.e. WITHOUT ANY WARRANTY of any kind, as
stated in the license.  In particular, its authors and contributors will take
no responsability for any possible bugs or any improper use of these programs,
including those resulting in incorrect scientific publications.
