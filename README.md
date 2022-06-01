TRIQS/SOM: Stochastic Optimization Method for Analytic Continuation
===================================================================

[![CI](https://github.com/krivenko/som/actions/workflows/CI.yml/badge.svg)](
https://github.com/krivenko/som/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-red)](
https://krivenko.github.io/som)

Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko @ gmail.com>

This program is an optimized implementation of an analytic continuation
method proposed by Andrey S. Mishchenko. A detailed description of
the method can be found in
<http://www.cond-mat.de/events/correl12/manuscripts/mishchenko.pdf>.

SOM 2.x offers new features, including the more advanced SOCC ([stochastic
optimization with consistent constraints](
https://doi.org/10.1103/PhysRevB.95.014102)) analytic continuation protocol.

SOM uses the TRIQS library version 3.1.x (SOM versions 1.x were based on
TRIQS 1.4.2). The TRIQS website is under <https://triqs.github.io/triqs/>.
Start there to learn about TRIQS.

A legacy documentation web-site for SOM 1.x is still available at
<https://krivenko.github.io/som1>.

Installation
------------

Please, refer to <http://krivenko.github.io/som/install.html> or
`doc/install.rst` in the source code tree for installation instructions.

Citing
------

[Accompanying Computer Physics Communications paper](
https://doi.org/10.1016/j.cpc.2019.01.021)
[[arXiv:1808.00603](https://arxiv.org/abs/1808.00603)]

```BibTeX
@article{SOM,
    title = {{TRIQS/SOM: Implementation of the stochastic optimization method
              for analytic continuation}},
    author = {Igor Krivenko and Malte Harland},
    journal = {Computer Physics Communications},
    volume = {239},
    pages = {166-183},
    year = {2019},
    issn = {0010-4655},
    doi = {https://doi.org/10.1016/j.cpc.2019.01.021},
    url = {https://www.sciencedirect.com/science/article/pii/S0010465519300402}
}
```

License
-------

SOM is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

SOM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
SOM (in the file LICENSE.txt in this directory). If not, see
<http://www.gnu.org/licenses/>.
