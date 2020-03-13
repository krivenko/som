SOM: Stochastic Optimization Method for Analytic Continuation
=============================================================

[![Build Status](https://travis-ci.org/krivenko/som.svg?branch=master)](https://travis-ci.org/krivenko/som)

http://krivenko.github.io/som/

[CPC paper](https://doi.org/10.1016/j.cpc.2019.01.021) [[arXiv:1808.00603](https://arxiv.org/abs/1808.00603)]

Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko @ gmail.com>

This program is an optimized implementation of an analytic continuation
method proposed by Andrey S. Mishchenko. A detailed description of
the method can be found under
http://www.cond-mat.de/events/correl12/manuscripts/mishchenko.pdf

This project uses the TRIQS library version 1.4.2.
The TRIQS website is under https://triqs.github.io/triqs/.
Start there to learn about TRIQS.

Installation
------------

Please, refer to http://krivenko.github.io/som/install.html or `doc/install.rst`
in the source code tree for installation instructions.

Source directory structure
--------------------------

 * `benchmark` - sets of Python scripts used to benchmark SOM's performance
   * `chi` - analytic continuation of QMC charge susceptibility
   * `microgap` - analytic continuation of a model spectrum with a very narrow gap
   * `semicircle` - analytic continuation of a semi-elliptic spectral function
   * `triangles` - analytic continuation of a spectrum made of two triangles
 * `c++` - main directory with C++ source files
   * `kernels` - C++ header files implementing individual integral kernels
   * `numerics` - C++ classes/routines for numerical algorithms (splines, Simpson's rule, dilogarithm function, etc)
 * `CHANGELOG.md` -  a brief list of changes between released versions
 * `cmake` - auxiliary files used by CMake to build application
 * `CMakeLists.txt` - main CMake configuration script
 * `LICENSE.txt` - text of the GNU General Public License
 * `python` - script files used by TRIQS Python wrapping tool to build the Python module of SOM
 * `test` - source code of unit tests executed by `make test`
   * `c++` - C++ unit tests
   * `CMakeLists.txt` - CMake configuration file for running unit tests
 * `doc` - documentation
   * `CMakeLists.txt` - CMake configuration file for building HTML documentation
   * `examples` - usage examples
   * `_static` - images used to build the HTML docs
   * `_templates` - auxiliary files used to build the HTML docs
   * `conf.py.in` - configuration file for Sphinx documentation generator
   * `notes` - LaTeX sources of implementation notes and Mathematica notebooks
 * `README.md` - this file
 * `.clang-format` - style file for clang-format utility
 * `.travis.yml` - configuration file for Travis CI continuous integration service

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
