##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# SOM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# SOM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# SOM. If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

from cpp2py.wrap_generator import module_

# The module
module = module_(full_name="spectral_stats",
                 doc=r"Statistical analysis of ensembles of spectral functions",
                 app_name="som")

# Imports
module.add_imports('som')

module.add_include("som/spectral_stats.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/vector.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
using namespace nda;
using namespace som;
""")

r_func_names = ["rectangle", "lorentzian", "gaussian"]
r_func_docstring = ', '.join([f"``{r_func}``" for r_func in r_func_names])

module.add_enum("resolution_function",
                r_func_names,
                "som",
                r"Type of resolution function :math:`\bar K(m, z)` used in spectral integrals.")

#
# spectral_integral()
#

module.add_function("double spectral_integral(double z_m, double delta_m, configuration c, resolution_function r_func)",
                    doc=r"""Evaluate spectral integral

.. math::

    i_m^{(j)} = \int_{-\infty}^\infty dz \bar K(m, z) A^{(j)}(z)

for a single energy interval.

**Parameters:**

:z_m: :class:`float`, Center of the energy interval.
:delta_m: :class:`float`, Length of the energy interval.
:c: :class:`som.Configuration`, Spectral function :math:`A^{(j)}(z)`.
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** :class:`float`, Value of the integral.
""" % r_func_docstring)

module.add_function("vector<double> spectral_integral(triqs::mesh::refreq mesh, configuration c, resolution_function r_func)",
                    doc=r"""Evaluate spectral integrals

.. math::

    i_m^{(j)} = \int_{-\infty}^\infty dz \bar K(m, z) A^{(j)}(z)

for energy intervals centered around points of a regular energy mesh.

**Parameters:**

:mesh: :class:`triqs.gf.meshes.MeshReFreq`, Real energy mesh.
:c: :class:`som.Configuration`, Spectral function :math:`A^{(j)}(z)`.
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** Real 1D NumPy array of values :math:`i_m^{(j)}`.
""" % r_func_docstring)

module.add_function("vector<double> spectral_integral(std::vector<std::pair<double, double>> intervals, configuration c, resolution_function r_func)",
                    doc=r"""Evaluate spectral integrals

.. math::

    i_m^{(j)} = \int_{-\infty}^\infty dz \bar K(m, z) A^{(j)}(z)

for a list of real energy intervals.

**Parameters:**

:intervals: :class:`list` [:class:`float`, :class:`float`], List of pairs
            (left interval boundary, right interval boundary).
:c: :class:`som.Configuration`, Spectral function :math:`A^{(j)}(z)`.
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** Real 1D NumPy array of values :math:`i_m^{(j)}`.
""" % r_func_docstring)

#
# spectral_avg()
#

module.add_function("vector<double> spectral_avg(som_core cont, long i, triqs::mesh::refreq mesh, resolution_function r_func)",
                    doc=r"""Compute spectral averages over a set of accumulated
particular solutions

.. math::

    i_m = \frac{1}{J} \sum_{j=1}^J i_m^{(j)}

for energy intervals centered around points of a regular energy mesh.

**Parameters:**

:cont: :class:`Som`, Analytic continuation object.
:i: :class:`int`, Index of the diagonal matrix element of the observable used
    to construct ``cont``.
:mesh: :class:`triqs.gf.meshes.MeshReFreq`, Real energy mesh.
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** Real 1D NumPy array of values :math:`i_m`.
""" % r_func_docstring)

module.add_function("vector<double> spectral_avg(som_core cont, long i, std::vector<std::pair<double, double>> intervals, resolution_function r_func)",
                    doc=r"""Compute spectral averages over a set of accumulated
particular solutions

.. math::

    i_m = \frac{1}{J} \sum_{j=1}^J i_m^{(j)}

for a list of real energy intervals.

**Parameters:**

:cont: :class:`Som`, Analytic continuation object.
:i: :class:`int`, Index of the diagonal matrix element of the observable used
    to construct ``cont``.
:intervals: :class:`list` [:class:`float`, :class:`float`], List of pairs
            (left interval boundary, right interval boundary).
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** Real 1D NumPy array of values :math:`i_m`.
""" % r_func_docstring)

#
# spectral_disp()
#

module.add_function("vector<double> spectral_disp(som_core cont, long i, triqs::mesh::refreq mesh, vector<double> avg, resolution_function r_func)",
                    doc=r"""Compute spectral dispersions of a set of accumulated
particular solutions

.. math::

    \sigma^2_m = \frac{1}{J} \sum_{j=1}^J (i_m^{(j)} - i_m)^2

for energy intervals centered around points of a regular energy mesh.

**Parameters:**

:cont: :class:`Som`, Analytic continuation object.
:i: :class:`int`, Index of the diagonal matrix element of the observable used
    to construct ``cont``.
:mesh: :class:`triqs.gf.meshes.MeshReFreq`, Real energy mesh.
:avg: Real 1D NumPy array of precomputed averages :math:`i_m`.
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** Real 1D NumPy array of values :math:`\sigma_m`.
""" % r_func_docstring)

module.add_function("vector<double> spectral_disp(som_core cont, long i, std::vector<std::pair<double, double>> intervals, vector<double> avg, resolution_function r_func)",
                    doc=r"""Compute spectral dispersions of a set of accumulated
particular solutions

.. math::

    \sigma^2_m = \frac{1}{J} \sum_{j=1}^J (i_m^{(j)} - i_m)^2

for a list of real energy intervals.

**Parameters:**

:cont: :class:`Som`, Analytic continuation object.
:i: :class:`int`, Index of the diagonal matrix element of the observable used
    to construct ``cont``.
:intervals: :class:`list` [:class:`float`, :class:`float`], List of pairs
            (left interval boundary, right interval boundary).
:avg: Real 1D NumPy array of precomputed averages :math:`i_m`.
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** Real 1D NumPy array of values :math:`\sigma_m`.
""" % r_func_docstring)

#
# spectral_corr()
#

module.add_function("matrix<double> spectral_corr(som_core cont, long i, triqs::mesh::refreq mesh, vector<double> avg, resolution_function r_func)",
                    doc=r"""Compute spectral two-point correlators of a set
of accumulated particular solutions

.. math::

    \sigma_{mm'} = \frac{1}{J} \sum_{j=1}^J (i_m^{(j)} - i_m)
                                            (i_{m'}^{(j)} - i_{m'})

for energy intervals centered around points of a regular energy mesh.

**Parameters:**

:cont: :class:`Som`, Analytic continuation object.
:i: :class:`int`, Index of the diagonal matrix element of the observable used
    to construct ``cont``.
:mesh: :class:`triqs.gf.meshes.MeshReFreq`, Real energy mesh.
:avg: Real 1D NumPy array of precomputed averages :math:`i_m`.
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** Real 2D NumPy array of values :math:`\sigma_{mm'}`.
""" % r_func_docstring)

module.add_function("matrix<double> spectral_corr(som_core cont, long i, std::vector<std::pair<double, double>> intervals, vector<double> avg, resolution_function r_func)",
                    doc=r"""Compute spectral two-point correlators of a set
of accumulated particular solutions

.. math::

    \sigma_{mm'} = \frac{1}{J} \sum_{j=1}^J (i_m^{(j)} - i_m)
                                            (i_{m'}^{(j)} - i_{m'})

for a list of real energy intervals.

**Parameters:**

:cont: :class:`Som`, Analytic continuation object.
:i: :class:`int`, Index of the diagonal matrix element of the observable used
    to construct ``cont``.
:intervals: :class:`list` [:class:`float`, :class:`float`], List of pairs
            (left interval boundary, right interval boundary).
:avg: Real 1D NumPy array of precomputed averages :math:`i_m`.
:r_func: :class:`str`, Name of the :ref:`resolution function
         <resolution_functions>` :math:`\bar K(m, z)`, one of %s.

**Returns:** Real 2D NumPy array of values :math:`\sigma_{mm'}`.
""" % r_func_docstring)

module.generate_code()
