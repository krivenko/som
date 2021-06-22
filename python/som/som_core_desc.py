##############################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
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
module = module_(full_name = "som_core",
                 doc = r"The Stochastic Optimization Method",
                 app_name = "som")

# Imports
module.add_imports('triqs.gf', 'triqs.statistics.histograms')

module.add_include("som/som_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/tuple.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
using namespace triqs::gfs;
using namespace triqs::statistics;
using namespace som;
""")

#
# Class Rectangle
#

c = class_(
    py_type = "Rectangle",
    c_type = "rectangle",
    is_printable = True,
    doc = r"Basis element of a spectral function"
)

c.add_property(name = "center",
               getter = cfunction(signature = "double()", calling_pattern = "auto result = self_c.center"),
               doc = "Center of the rectangle")
c.add_property(name = "width",
               getter = cfunction(signature = "double()", calling_pattern = "auto result = self_c.width"),
               doc = "Width of the rectangle")
c.add_property(name = "height",
               getter = cfunction(signature = "double()", calling_pattern = "auto result = self_c.height"),
               doc = "Height of the rectangle")

c.add_property(name = "norm",
               getter = cfunction("double norm()"),
               doc = "Norm (area) of the rectangle")

c.add_call(signature = "double(double x)",
           doc = "Substitute value 'x' into the rectangle function")

module.add_class(c)

#
# Class Configuration
#

c = class_(
    py_type = "Configuration",
    c_type = "configuration",
    is_printable = True,
    hdf5 = True,
    doc = r"Sum of rectangles"
)

c.add_len(doc = "Number of rectangles in the configuration")
c.add_getitem(signature = "rectangle(int i)",
              calling_pattern = """
                  int len = self_c.size();
                  if((i < -len) || (i >= len)) CPP2PY_RUNTIME_ERROR << "Rectangle index " << i << " is out of bounds";
                  auto const& result = self_c[i >= 0 ? i : len + i];
              """,
              doc = "Individual rectangle access")
c.add_iterator()

c.add_property(name = "norm",
               getter = cfunction("double norm()"),
               doc = "Norm of the configuration")

module.add_class(c)

#
# Class SomCore
#

module.add_enum("observable_kind", ["FermionGf","BosonCorr","BosonAutoCorr","ZeroTemp"], "som", "Kinds of observables")

c = class_(
    py_type = "SomCore",        # name of the python class
    c_type = "som_core",        # name of the C++ class
    doc = r"Main class of SOM"  # doc of the C++ class
)

c.add_constructor("""(gf_view<imtime> g_tau, gf_view<imtime> S_tau, som::observable_kind kind = FermionGf, vector<double> norms = {})""",
                  doc = """Construct on imaginary-time quantities """)

c.add_constructor("""(gf_view<imfreq> g_iw, gf_view<imfreq> S_iw, som::observable_kind kind = FermionGf, vector<double> norms = {})""",
                  doc = """Construct on imaginary-frequency quantities """)

c.add_constructor("""(gf_view<legendre> g_l, gf_view<legendre> S_l, som::observable_kind kind = FermionGf, vector<double> norms = {})""",
                  doc = """Construct on quantities in Legendre polynomial basis """)

# TODO: Update docstrings

c.add_method("""void run (**som::run_parameters_t)""",
             release_GIL_and_enable_signal = True,
             doc = """
Main parameters
---------------

+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+
| Parameter Name  | Type          | Default                       | Documentation                                                                                            |
+=================+===============+===============================+==========================================================================================================+
| energy_window   | (float,float) | --                            | Estimated lower and upper bounds of the spectrum.                                                        |
|                 |               |                               | Negative values of the lower bound will be reset to 0 for susceptibilities and conductivity.             |
+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+
| max_time        | int           | -1 = infinite                 | Maximum runtime in seconds, use -1 to set infinite.                                                      |
+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+
| verbosity       | int           | 2 on MPI rank 0, 0 otherwise. | Verbosity level (max level - 3).                                                                         |
+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+
| t               | int           | 50                            | Number of elementary updates per global update (:math:`T`).                                              |
|                 |               |                               | Bigger values make the algorithm more ergodic.                                                           |
+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+
| adjust_f        | bool          | False                         | Adjust the number of global updates automatically.                                                       |
+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+
| l               | int           | 2000                          | Number of particular solutions used in the final accumulation (:math:`L`); ignored if `adjust_l = True`. |
|                 |               |                               | Bigger values reduce noise in the final solution / make it smoother.                                     |
+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+
| adjust_l        | bool          | False                         | Adjust the number of solutions used in the final accumulation.                                           |
+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+
| make_histograms | bool          | False                         | Accumulate histograms of objective function values.                                                      |
+-----------------+---------------+-------------------------------+----------------------------------------------------------------------------------------------------------+

Fine tuning options
-------------------

+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| Parameter Name      | Type      | Default                       | Documentation                                                                                       |
+=====================+===========+===============================+=====================================================================================================+
| random_seed         | int       | 34788 + 928374 * MPI.rank     | Seed for random number generator.                                                                   |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| random_name         | str       | ""                            | Name of random number generator.                                                                    |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| max_rects           | int       | 60                            | Maximum number of rectangles to represent spectra (:math:`K_{max}`), should be below 70.            |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| min_rect_width      | float     | 1e-3                          | Minimal width of a rectangle, in units of the energy window width.                                  |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| min_rect_weight     | float     | 1e-3                          | Minimal weight of a rectangle, in units of the requested solution norm.                             |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| distrib_d_max       | float     | 2                             | Maximal parameter of the power-law distribution function for the Metropolis algorithm.              |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| gamma               | float     | 2                             | Proposal probability parameter :math:`\gamma`.                                                      |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| adjust_f_range      | (int,int) | (100,5000)                    | Search range for the number of global updates.                                                      |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| adjust_f_l          | int       | 20                            | Number of particular solutions used to adjust :math:`F`.                                            |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| adjust_f_kappa      | float     | 0.25                          | Limiting value of :math:`\kappa` used to adjust :math:`F`.                                          |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| adjust_l_range      | (int,int) | (100,2000)                    | Search range for the number of solutions used in the final accumulation.                            |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| adjust_l_good_d     | float     | 2.0                           | Maximal ratio :math:`D/D_\mathrm{min}` for a particular solution to be considered good.             |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| adjust_l_verygood_d | float     | 4/3                           | Maximal ratio :math:`D/D_\mathrm{min}` for a particular solution to be considered very good.        |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| adjust_l_ratio      | float     | 0.95                          | Critical ratio :math:`N_\mathrm{very good}/N_\mathrm{good}` to stop :math:`L`-adjustment procedure. |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| hist_max            | float     | 2.0                           | Right boundary of the histograms, in units of :math:`D_\mathrm{min}`                                |
|                     |           |                               | (left boundary is always set to :math:`D_\mathrm{min}`).                                            |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
| hist_n_bins         | int       | 100                           | Number of bins for the histograms.                                                                  |
+---------------------+-----------+-------------------------------+-----------------------------------------------------------------------------------------------------+
""")

#
# Add parameters from worker_parameters_t
#

def add_worker_parameters(conv):
    conv.add_member(c_name = "energy_window",
                    c_type = "std::pair<double,double>",
                    initializer = """  """,
                    doc = """Estimated lower and upper bounds of the spectrum.\nNegative values of the lower bound will be reset to 0 for susceptibilities and conductivity.""")

    conv.add_member(c_name = "max_time",
                    c_type = "int",
                    initializer = """ -1 """,
                    doc = """Maximum runtime in seconds, use -1 to set infinite.""")

    conv.add_member(c_name = "verbosity",
                    c_type = "int",
                    initializer = """ ((mpi::communicator().rank() == 0) ? 2 : 0) """,
                    doc = """Verbosity level (max level - 3).""")

    conv.add_member(c_name = "t",
                    c_type = "int",
                    initializer = """50""",
                    doc = """Number of elementary updates per global update (:math:`T`).\nBigger values make the algorithm more ergodic.""")

    conv.add_member(c_name = "random_seed",
                    c_type = "int",
                    initializer = """34788 + 928374 * mpi::communicator().rank()""",
                    doc = """Seed for random number generator.""")

    conv.add_member(c_name = "random_name",
                    c_type = "std::string",
                    initializer = """ "" """,
                    doc = """Name of random number generator.""")

    conv.add_member(c_name = "max_rects",
                    c_type = "int",
                    initializer = """60""",
                    doc = """Maximum number of rectangles to represent spectra (:math:`K_{max}`), should be below 70.""")

    conv.add_member(c_name = "min_rect_width",
                    c_type = "double",
                    initializer = """1e-3""",
                    doc = """Minimal width of a rectangle, in units of the energy window width.""")

    conv.add_member(c_name = "min_rect_weight",
                    c_type = "double",
                    initializer = """1e-3""",
                    doc = """Minimal weight of a rectangle, in units of the requested solution norm.""")

    conv.add_member(c_name = "distrib_d_max",
                    c_type = "double",
                    initializer = """2""",
                    doc = """Maximal parameter of the power-law distribution function for the Metropolis algorithm.""")

    conv.add_member(c_name = "gamma",
                    c_type = "double",
                    initializer = """2""",
                    doc = """Proposal probability parameter :math:`\gamma`.""")

#
# Converter for run_parameters_t
#

run_params_conv = converter_(
    c_type = "som::run_parameters_t",
    doc = """Arguments of SomCore.run()""",
)

add_worker_parameters(run_params_conv)

run_params_conv.add_member(c_name = "t",
                           c_type = "int",
                           initializer = """50""",
                           doc = """Number of elementary updates per global update (:math:`T`).\nBigger values make the algorithm more ergodic.""")

run_params_conv.add_member(c_name = "f",
                           c_type = "int",
                           initializer = """100""",
                           doc = """Number of global updates (:math:`F`); ignored if `adjust_f = True`.\nBigger values make the algorithm more ergodic.""")

run_params_conv.add_member(c_name = "l",
                           c_type = "int",
                           initializer = """2000""",
                           doc = """Number of particular solutions used in the final accumulation (:math:`L`); ignored if `adjust_l = True`.\nBigger values reduce noise in the final solution / make it smoother.""")

run_params_conv.add_member(c_name = "adjust_l",
                           c_type = "bool",
                           initializer = """false""",
                           doc = """Adjust the number of solutions used in the final accumulation.""")

run_params_conv.add_member(c_name = "make_histograms",
                           c_type = "bool",
                           initializer = """false""",
                           doc = """Accumulate histograms of objective function values.""")

run_params_conv.add_member(c_name = "adjust_l_range",
                           c_type = "std::pair<int,int>",
                           initializer = """std::pair<int,int>{100,2000}""",
                           doc = """Search range for the number of solutions used in the final accumulation.""")

run_params_conv.add_member(c_name = "adjust_l_good_d",
                           c_type = "double",
                           initializer = """2.0""",
                           doc = """Maximal ratio :math:`D/D_\mathrm{min}` for a particular solution to be considered good.""")

run_params_conv.add_member(c_name = "adjust_l_verygood_d",
                           c_type = "double",
                           initializer = """4.0/3.0""",
                           doc = """Maximal ratio :math:`D/D_\mathrm{min}` for a particular solution to be considered very good.""")

run_params_conv.add_member(c_name = "adjust_l_ratio",
                           c_type = "double",
                           initializer = """0.95""",
                           doc = """Critical ratio :math:`N_\mathrm{very good}/N_\mathrm{good}` to stop :math:`L`-adjustment procedure.""")

run_params_conv.add_member(c_name = "hist_max",
                           c_type = "double",
                           initializer = """2.0""",
                           doc = """Right boundary of the histograms, in units of :math:`D_\mathrm{min}`\n(left boundary is always set to :math:`D_\mathrm{min}`).""")

run_params_conv.add_member(c_name = "hist_n_bins",
                           c_type = "int",
                           initializer = """100""",
                           doc = """Number of bins for the histograms.""")

module.add_converter(run_params_conv)

#
# SomCore.fill_observable()
#

c.add_method("void fill_observable(gf_view<refreq> g_w)", calling_pattern = "g_w() = self_c")
c.add_method("void fill_observable(gf_view<imtime> g_tau)", calling_pattern = "g_tau() = self_c")
c.add_method("void fill_observable(gf_view<imfreq> g_iw)", calling_pattern = "g_iw() = self_c")
c.add_method("void fill_observable(gf_view<legendre> g_l)", calling_pattern = "g_l() = self_c")

#
# SomCore.solution() and SomCore.solutions
#

c.add_method(name = "solution",
             signature = "configuration get_solution (int i)",
             doc = """Final solution for the i-th diagonal matrix element of the observable""")

c.add_property(name = "solutions",
               getter = cfunction("std::vector<configuration> get_solutions ()"),
               doc = """Final solutions, one per diagonal matrix element of the observable""")

#
# SomCore.histogram() and SomCore.histograms
#

c.add_method(name = "histogram",
             signature = "std::optional<histogram> get_histogram(int i)",
             doc = """Accumulated objective function histogram for the i-th diagonal matrix element of the observable""")

c.add_property(name = "histograms",
               getter = cfunction(signature = "std::optional<std::vector<histogram>> get_histograms ()"),
               doc = """Accumulated objective function histograms, one per diagonal matrix element of the observable""")

#
# SomCore.compute_tail()
#

c.add_method("triqs::arrays::array<dcomplex, 3> compute_tail(int max_order)",
             doc = """Compute GF tail coefficients using calculated spectra""")

#
# Other attributes of SomCore
#

c.add_property(name = "last_run_parameters",
               getter = cfunction("som::run_parameters_t get_last_run_parameters ()"),
               doc = """Set of parameters used in the last call to run() """)

c.add_property(name = "run_status",
               getter = cfunction("int get_run_status ()"),
               doc = """Status of the run on exit """)

#
# Converter for adjust_f_parameters_t
#

adjust_f_params_conv = converter_(
    c_type = "som::adjust_f_parameters_t",
    doc = """Arguments of adjust_f()""",
)

add_worker_parameters(adjust_f_params_conv)

adjust_f_params_conv.add_member(c_name = "f_range",
                                c_type = "std::pair<int,int>",
                                initializer = """std::pair<int,int>{100,5000}""",
                                doc = """Search range for the number of global updates.""")

adjust_f_params_conv.add_member(c_name = "l",
                                c_type = "int",
                                initializer = """20""",
                                doc = """Number of particular solutions used to adjust :math:`F`.""")

adjust_f_params_conv.add_member(c_name = "kappa",
                                c_type = "double",
                                initializer = """0.25""",
                                doc = """Limiting value of :math:`\kappa` used to adjust :math:`F`.""")

module.add_converter(adjust_f_params_conv)

#
# SomCore.adjust_f()
#

c.add_method("int adjust_f(**som::adjust_f_parameters_t)",
             doc = """Automatically adjust the number of global updates""")

module.add_class(c)

module.generate_code()
