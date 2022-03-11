##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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

from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "spectral_stats",
                 doc = r"Statistical analysis of noisy spectral functions",
                 app_name = "som")

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

module.add_enum("resolution_function",
                ["rectangle", "lorentzian", "gaussian"],
                "som",
                r"Type of resolution function :math:`\bar K(m, z)` used in spectral integrals.")

#
# spectral_integral()
#

module.add_function("double spectral_integral(double z_m, double delta_m, configuration c, resolution_function r_func)",
                    doc = r"""Spectral integral :math:`i_m^{(j)} = \int_{-\infty}^\infty dz \bar K(m, z) A^{(j)}(z)`.""")

module.add_function("vector<double> spectral_integral(triqs::mesh::refreq mesh, configuration c, resolution_function r_func)",
                    doc = r"""Spectral integrals :math:`i_m^{(j)} = \int_{-\infty}^\infty dz \bar K(m, z) A^{(j)}(z)`\
parametrized by real energy points from a regular mesh.""")

module.add_function("vector<double> spectral_integral(std::vector<std::pair<double, double>> intervals, configuration c, resolution_function r_func)",
                    doc = r"""Spectral integrals :math:`i_m^{(j)} = \int_{-\infty}^\infty dz \bar K(m, z) A^{(j)}(z)`
parametrized by a list of real energy intervals.""")

#
# spectral_avg()
#

module.add_function("vector<double> spectral_avg(som_core cont, int i, triqs::mesh::refreq mesh, resolution_function r_func)",
                    doc = r"""Compute spectral averages :math:`i_m = J^{-1} \sum_{j=1}^J i_m^{(j)}`
parametrized by real energy points from a regular mesh.""")

module.add_function("vector<double> spectral_avg(som_core cont, int i, std::vector<std::pair<double, double>> intervals, resolution_function r_func)",
                    doc = r"""Compute spectral averages :math:`i_m = J^{-1} \sum_{j=1}^J i_m^{(j)}`
parametrized by a list of real energy intervals.""")

#
# spectral_disp()
#

module.add_function("vector<double> spectral_disp(som_core cont, int i, triqs::mesh::refreq mesh, vector<double> avg, resolution_function r_func)",
                    doc = r"""Compute spectral dispersions :math:`\sigma^2_m = J^{-1} \sum_{j=1}^J (i_m^{(j)} - i_m)^2`
parametrized by real energy points from a regular mesh.""")

module.add_function("vector<double> spectral_disp(som_core cont, int i, std::vector<std::pair<double, double>> intervals, vector<double> avg, resolution_function r_func)",
                    doc = r"""Compute spectral dispersions :math:`\sigma^2_m = J^{-1} \sum_{j=1}^J (i_m^{(j)} - i_m)^2`
parametrized by a list of real energy intervals.""")

#
# spectral_corr()
#

module.add_function("matrix<double> spectral_corr(som_core cont, int i, triqs::mesh::refreq mesh, vector<double> avg, resolution_function r_func)",
                    doc = r"""Compute spectral correlators :math:`\sigma_{mm'} = J^{-1} \sum_{j=1}^J (i_m^{(j)} - i_m) (i_{m'}^{(j)} - i_{m'})`
parametrized by real energy points from a regular mesh.""")

module.add_function("matrix<double> spectral_corr(som_core cont, int i, std::vector<std::pair<double, double>> intervals, vector<double> avg, resolution_function r_func)",
                    doc = r"""Compute spectral correlators :math:`\sigma_{mm'} = J^{-1} \sum_{j=1}^J (i_m^{(j)} - i_m) (i_{m'}^{(j)} - i_{m'})`
parametrized by a list of real energy intervals.""")

module.generate_code()
