TRIQS/SOM: Stochastic Optimization Method for Analytic Continuation
===================================================================

[![CI](https://github.com/krivenko/som/actions/workflows/build-and-deploy.yml/badge.svg)](
https://github.com/krivenko/som/actions/workflows/build-and-deploy.yml)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-red)](
https://krivenko.github.io/som)

Copyright (C) 2016-2026 Igor Krivenko <iskrivenko [at] proton [dot] me>

This program is an optimized implementation of an analytic continuation
method proposed by Andrey S. Mishchenko. A detailed description of
the method can be found in
<http://www.cond-mat.de/events/correl12/manuscripts/mishchenko.pdf>.

SOM 2.x offers new features, including the more advanced SOCC ([stochastic
optimization with consistent constraints](
https://doi.org/10.1103/PhysRevB.95.014102)) analytic continuation protocol.

SOM uses the TRIQS library version 3.2.x or 3.3.x. The TRIQS website is under
<https://triqs.github.io/triqs/>. Start there to learn about TRIQS.

**Legacy SOM versions 1.x are based on TRIQS 1.4.2. The 1.x source code is still
available on the `1.x` branch; Its respective documentation web-site
was moved to <https://krivenko.github.io/som1>**.

Installation
------------

SOM can be installed using the [`conda`](https://www.anaconda.com/) package
manager or built from source either manually or by means of the
[EasyBuild](https://easybuild.io/) HPC software management framework.

We also offer prebuilt [Docker images](https://hub.docker.com/r/ikrivenko/som/).

Please, refer to <http://krivenko.github.io/som/install.html> or
`doc/install.rst` in the source code tree for the detailed installation
instructions.

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

[Version 2.0 announcement](https://doi.org/10.1016/j.cpc.2022.108491)

```BibTeX
@article{SOM2,
  title = {{TRIQS/SOM 2.0: Implementation of the stochastic optimization
            with consistent constraints for analytic continuation}},
  author = {Igor Krivenko and Andrey S. Mishchenko},
  journal = {Computer Physics Communications},
  volume = {280},
  pages = {108491},
  year = {2022},
  issn = {0010-4655},
  doi = {https://doi.org/10.1016/j.cpc.2022.108491},
  url = {https://www.sciencedirect.com/science/article/pii/S0010465522002107}
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
