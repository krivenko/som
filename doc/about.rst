.. _about:

About SOM
=========

.. highlight:: none

Authors
-------

This analytic continuation package is written by Igor Krivenko
(companion paper [KH2019]_,
`arXiv:1808.00603 <https://arxiv.org/abs/1808.00603>`_).

It is based on a method devised by Andrey S. Mishchenko and coauthors
[MPSS2000]_. You can find a detailed description of the method in Chapter 14 of
the Lecture Notes from the Autumn School "Correlated Electrons: From Models to
Materials".

::

    Eva Pavarini, Erik Koch, Frithjof Anders, and Mark Jarrell (eds.)
    Correlated Electrons: From Models to Materials
    Modeling and Simulation, Vol. 2
    Verlag des Forschungszentrum Jülich, 2012
    ISBN 978-3-89336-796-2

The Lecture Notes are available free of charge online,

    http://www.cond-mat.de/events/correl12/manuscripts/

The Stochastic Optimization with Consistent Constraints extensions are proposed
by Olga Goulko and coauthors in [GMPPS2017]_.

.. [KH2019] `I. Krivenko, M. Harland,
   Comput. Phys. Commun. 239, 166-183 (2019)
   <https://doi.org/10.1016/j.cpc.2019.01.021>`_
   (:download:`bibtex file <TRIQSSOM.bib>`)

.. [MPSS2000]
   `A. S. Mishchenko, N. V. Prokof'ev, A. Sakamoto and B. V. Svistunov,
   Phys. Rev. B 62, 6317 (2000) <https://doi.org/10.1103/PhysRevB.62.6317>`_
   (:download:`bibtex file <SOM.bib>`)

.. [GMPPS2017]
   `O. Goulko, A. S. Mishchenko, L. Pollet, N. Prokof'ev and B. Svistunov,
   Phys. Rev. B 95, 014102 (2017) <https://doi.org/10.1103/PhysRevB.95.014102>`_
   (:download:`bibtex file <SOCC.bib>`)

TRIQS
-----

It is highly recommended to get familiar with basic usage of the TRIQS library,
before starting to work with this program.

https://triqs.github.io/

Here is a simple tutorial showing how to manipulate Green's function objects of
TRIQS in Python scripts,

https://triqs.github.io/triqs/latest/userguide/gfs/gfs_tutorial_python.html

Acknowledgements
----------------

.. sidebar::

      .. image:: _static/logo_uhh.svg
         :width: 75%
         :align: center
         :target: https://www.uni-hamburg.de/en.html
         :alt: Universität Hamburg

      |

      .. image:: _static/logo_sfb668.jpg
         :width: 75%
         :align: center
         :target: http://www.sfb668.de/
         :alt: Sonderforschungsbereich 668

I am very grateful to Malte Harland for extensive and thorough testing of the
code at late development stages.

I would like to acknowledge support from the University of Hamburg, as well as
from the Deutsche Forschungsgemeinschaft via the Project SFB-668
"Magnetismus vom Einzelatom zur Nanostruktur" (A3). Some large-scale tests have
been run on the JURECA HPC machine of the Forschungszentrum Jülich (project
"Continuous Time Quantum Monte Carlo for materials").

Credits to Snir Gazit, who suggested adding support for full covariance matrices
in SOM 2.0.0.

License
-------

The SOM package is published under the `GNU General Public License, version 3
<http://www.gnu.org/licenses/gpl.html>`_.

Note that it *implies* that applications using SOM must also be GPL.

Usage disclaimer
----------------

The program is provided as is, i.e. WITHOUT ANY WARRANTY of any kind, as
stated in the license. In particular, its author and contributors will take
no responsibility for any possible bugs or any improper use of these programs,
including those resulting in incorrect scientific publications.
