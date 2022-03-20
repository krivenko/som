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
module = module_(full_name = "som_core",
                 doc = r"The Stochastic Optimization Method",
                 app_name = "som")

# Imports
module.add_imports('triqs.gf', 'triqs.stat.histograms')

module.add_include("som/som_core/som_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/tuple.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
using namespace triqs::gfs;
using namespace triqs::stat;
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
# Common fragments of docstrings
#

docstring_params_header_main = """
Main parameters
---------------
""".strip()

docstring_params_header_fine = """
Fine tuning options
-------------------
""".strip()

docstring_params_table_header = """
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| Parameter Name        | Type          | Default                       | Documentation                                                                                         |
+=======================+===============+===============================+=======================================================================================================+
""".strip()

docstring_worker_params_main = """
| energy_window         | (float,float) | --                            | Estimated lower and upper bounds of the spectrum.                                                     |
|                       |               |                               | Negative values of the lower bound will be reset to 0 for BosonAutoCorr and ZeroTemp observables.     |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| max_time              | int           | -1 = infinite                 | Maximum runtime in seconds, use -1 to set infinite.                                                   |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| verbosity             | int           | 2 on MPI rank 0, 0 otherwise. | Verbosity level (max level - 3).                                                                      |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| t                     | int           | 50                            | Number of elementary updates per global update (:math:`T`).                                           |
|                       |               |                               | Bigger values make the algorithm more ergodic.                                                        |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
""".strip()

docstring_worker_params_fine = """
| random_seed           | int           | 34788 + 928374 * MPI.rank     | Seed for random number generator.                                                                     |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| random_name           | str           | ""                            | Name of random number generator (MT19937 by default).                                                 |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| max_rects             | int           | 60                            | Maximum number of rectangles in a particular solution (:math:`K_{max}`), should be below 70.          |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| min_rect_width        | float         | 1e-3                          | Minimal width of a rectangle, in units of the energy window width.                                    |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| min_rect_weight       | float         | 1e-3                          | Minimal weight of a rectangle, in units of the requested solution norm.                               |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| distrib_d_max         | float         | 2                             | Maximal parameter of the power-law distribution function for the Metropolis algorithm.                |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| gamma                 | float         | 2                             | Proposal probability parameter :math:`\gamma`.                                                        |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
""".strip()

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
             doc = fr"""
Automatically adjust the number of global updates :math:`F`
===========================================================

{docstring_params_header_main}

{docstring_params_table_header}
{docstring_worker_params_main}
| f_range               | (int,int)     | (100,5000)                    | Search range for the number of global updates.                                                        |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+

{docstring_params_header_fine}

{docstring_params_table_header}
{docstring_worker_params_fine}
| l                     | int           | 20                            | Number of particular solutions used to adjust :math:`F`.                                              |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| kappa                 | float         | 0.25                          | Limiting value of :math:`\kappa` used to adjust :math:`F`.                                            |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
""")

#
# Converter for accumulate_parameters_t
#

accumulate_params_conv = converter_(
    c_type = "som::accumulate_parameters_t",
    doc = """Arguments of SomCore.accumulate()""",
)

add_worker_parameters(accumulate_params_conv)

accumulate_params_conv.add_member(c_name = "t",
                                  c_type = "int",
                                  initializer = """50""",
                                  doc = """Number of elementary updates per global update (:math:`T`).\nBigger values make the algorithm more ergodic.""")

accumulate_params_conv.add_member(c_name = "f",
                                  c_type = "int",
                                  initializer = """100""",
                                  doc = """Number of global updates (:math:`F`); ignored if `adjust_f = True`.\nBigger values make the algorithm more ergodic.""")

accumulate_params_conv.add_member(c_name = "l",
                                  c_type = "int",
                                  initializer = """2000""",
                                  doc = """Number of particular solutions used in the final accumulation (:math:`L`); ignored if `adjust_l = True`.\nBigger values reduce noise in the final solution / make it smoother.""")

accumulate_params_conv.add_member(c_name = "adjust_l",
                                  c_type = "bool",
                                  initializer = """false""",
                                  doc = """Adjust the number of solutions used in the final accumulation.""")

accumulate_params_conv.add_member(c_name = "make_histograms",
                                  c_type = "bool",
                                  initializer = """false""",
                                  doc = """Accumulate histograms of :math:`\chi` values.""")

accumulate_params_conv.add_member(c_name = "adjust_l_range",
                                  c_type = "std::pair<int,int>",
                                  initializer = """std::pair<int,int>{100,2000}""",
                                  doc = """Search range for the number of solutions used in the final accumulation.""")

accumulate_params_conv.add_member(c_name = "adjust_l_good_chi",
                                  c_type = "double",
                                  initializer = """2.0""",
                                  doc = r"""Maximal ratio :math:`\chi/\chi_\mathrm{min}` for a particular solution to be considered good.""")

accumulate_params_conv.add_member(c_name = "adjust_l_verygood_chi",
                                  c_type = "double",
                                  initializer = """4.0/3.0""",
                                  doc = r"""Maximal ratio :math:`\chi/\chi_\mathrm{min}` for a particular solution to be considered very good.""")

accumulate_params_conv.add_member(c_name = "adjust_l_ratio",
                                  c_type = "double",
                                  initializer = """0.95""",
                                  doc = r"""Critical ratio :math:`N_\mathrm{very good}/N_\mathrm{good}` to stop :math:`L`-adjustment procedure.""")

accumulate_params_conv.add_member(c_name = "hist_max",
                                  c_type = "double",
                                  initializer = """2.0""",
                                  doc = r"""Right boundary of the histograms, in units of :math:`\chi_\mathrm{min}`\n(left boundary is always set to :math:`\chi_\mathrm{min}`).""")

accumulate_params_conv.add_member(c_name = "hist_n_bins",
                                  c_type = "int",
                                  initializer = """100""",
                                  doc = """Number of bins for the histograms.""")

module.add_converter(accumulate_params_conv)

#
# SomCore.accumulate()
#

docstring_accumulate = fr"""
{docstring_params_header_main}

{docstring_params_table_header}
{docstring_worker_params_main}
| l                     | int           | 2000                          | Number of particular solutions to accumulate (:math:`L`); ignored if `adjust_l = True`.               |
|                       |               |                               | Bigger values reduce noise in the final solution / make it smoother.                                  |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| adjust_l              | bool          | False                         | Automatically adjust the number of particular solutions to accumulate.                                |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| make_histograms       | bool          | False                         | Accumulate histograms of :math:`\chi` values.                                                         |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+

{docstring_params_header_fine}

{docstring_params_table_header}
{docstring_worker_params_fine}
| adjust_l_range        | (int,int)     | (100,2000)                    | Search range for the number of particular solutions to accumulate.                                    |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| adjust_l_good_chi     | float         | 2.0                           | Maximal ratio :math:`\chi/\chi_\mathrm{{min}}` for a particular solution to be considered good.         |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| adjust_l_verygood_chi | float         | 4/3                           | Maximal ratio :math:`\chi/\chi_\mathrm{{min}}` for a particular solution to be considered very good.    |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| adjust_l_ratio        | float         | 0.95                          | Critical ratio :math:`N_\mathrm{{very good}}/N_\mathrm{{good}}` to stop                                   |
|                       |               |                               | the :math:`L`-adjustment procedure.                                                                   |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| hist_max              | float         | 2.0                           | Right boundary of the histograms, in units of :math:`\chi_\mathrm{{min}}`                               |
|                       |               |                               | (left boundary is always set to :math:`\chi_\mathrm{{min}}`).                                           |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| hist_n_bins           | int           | 100                           | Number of bins for the histograms.                                                                    |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
"""

c.add_method("""void accumulate(**som::accumulate_parameters_t)""",
             release_GIL_and_enable_signal = True,
             doc = f"""
Accumulate particular solutions
===============================

{docstring_accumulate}
""")

#
# SomCore.compute_final_solution()
#

c.add_method("std::vector<double> compute_final_solution(double good_chi_rel = 2.0, double good_chi_abs = HUGE_VAL)",
             doc = """Select particular solutions according to the standard SOM criterion and compute the final solution""")

#
# Converter for final_solution_cc_parameters_t
#

compute_final_solution_cc_params_conv = converter_(
    c_type = "som::final_solution_cc_parameters_t",
    doc = """Arguments of SomCore.compute_final_solution_cc()""",
)

compute_final_solution_cc_params_conv.add_member(c_name = "refreq_mesh",
                                                 c_type = "std::variant<triqs::mesh::refreq, nda::array<double, 1>>",
                                                 initializer = """ """,
                                                 doc = """Grid of energy points used in derivative regularization procedure.""")

compute_final_solution_cc_params_conv.add_member(c_name = "verbosity",
                                                 c_type = "int",
                                                 initializer = """ ((mpi::communicator().rank() == 0) ? 1 : 0) """,
                                                 doc = """Verbosity level (max level - 2).""")

compute_final_solution_cc_params_conv.add_member(c_name = "good_chi_rel",
                                                 c_type = "double",
                                                 initializer = "2.0",
                                                 doc = """Maximal ratio :math:`\chi/\chi_\mathrm{min}` for a particular solution to be selected.
This criterion must be fulfilled together with the one set by `good_chi_abs`.""")

compute_final_solution_cc_params_conv.add_member(c_name = "good_chi_abs",
                                                 c_type = "double",
                                                 initializer = "HUGE_VAL",
                                                 doc = """Maximal value of :math:`\chi` for a particular solution to be selected.
This criterion must be fulfilled together with the one set by `good_chi_rel`.""")

compute_final_solution_cc_params_conv.add_member(c_name = "default_model",
                                                 c_type = "nda::array<double, 1>",
                                                 initializer = """{}""",
                                                 doc = """Default model of the spectral function evaluated at energy points of `refreq_mesh`.""")

compute_final_solution_cc_params_conv.add_member(c_name = "default_model_weights",
                                                 c_type = "nda::array<double, 1>",
                                                 initializer = """{}""",
                                                 doc = """Weights determining how much deviations from `default_model` are penalized at each energy point of `refreq_mesh`.""")

compute_final_solution_cc_params_conv.add_member(c_name = "max_iter",
                                                 c_type = "int",
                                                 initializer = """20""",
                                                 doc = """Maximum allowed number of parameter adjustment iterations.""")
compute_final_solution_cc_params_conv.add_member(c_name = "unity_sum_coeff",
                                                 c_type = "double",
                                                 initializer = """1e6""",
                                                 doc = """Coefficient of the term that enforces the unity sum constraint.""")

compute_final_solution_cc_params_conv.add_member(c_name = "amp_penalty_max",
                                                 c_type = "double",
                                                 initializer = """1e6""",
                                                 doc = """Maximum value of the regularization parameter that penalizes negative values of the spectral function.""")

compute_final_solution_cc_params_conv.add_member(c_name = "amp_penalty_divisor",
                                                 c_type = "double",
                                                 initializer = """100""",
                                                 doc = """Divisor used to reduce the regularization parameter that penalizes negative values of the spectral function.""")

compute_final_solution_cc_params_conv.add_member(c_name = "der_penalty_init",
                                                 c_type = "double",
                                                 initializer = """1.0""",
                                                 doc = """Initial value of the regularization parameters that penalize large derivatives of the solution.""")

compute_final_solution_cc_params_conv.add_member(c_name = "der_penalty_coeff",
                                                 c_type = "double",
                                                 initializer = """2.0""",
                                                 doc = """Coefficient used to increase the regularization parameters that penalize large derivatives of the solution.""")

module.add_converter(compute_final_solution_cc_params_conv)

#
# SomCore.compute_final_solution_cc()
#

c.add_method("std::vector<double> compute_final_solution_cc(**som::final_solution_cc_parameters_t)",
             #release_GIL_and_enable_signal = True, FIXME
             doc = f"""
Compute the final solution using the SOCC protocol
==================================================

{docstring_params_header_main}

{docstring_params_table_header}
| refreq_mesh           | MeshReFreq or | --                            | Grid of energy points :math:`\epsilon_k` used in CC regularization procedure.                         |
|                       | real 1D       |                               | Either a TRIQS real-frequency mesh object or a strictly ordered list of points is allowed.            |
|                       | array_like    |                               |                                                                                                       |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| verbosity             | int           | 1 on MPI rank 0, 0 otherwise. | Verbosity level (max level - 2).                                                                      |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| good_chi_rel          | float         | 2.0                           | Maximal ratio :math:`\chi/\chi_\mathrm{{min}}` for a particular solution to be selected.                |
|                       |               |                               | This criterion must be fulfilled together with the one set by `good_chi_abs`.                         |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| good_chi_abs          | float         | infinity                      | Maximal value of :math:`\chi` for a particular solution to be selected.                               |
|                       |               |                               | This criterion must be fulfilled together with the one set by `good_chi_rel`.                         |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| default_model         | Real 1D array | []                            | Optional default model of the spectral function evaluated at energy points of `refreq_mesh`           |
|                       |               |                               | (:math:`A_T(\epsilon_k)`).                                                                            |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| default_model_weights | Real 1D array | []                            | Weights determining how much deviations from `default_model` are penalized at each energy point       |
|                       |               |                               | (:math:`T_k`).                                                                                        |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| max_iter              | int           | 20                            | Maximum allowed number of parameter adjustment iterations.                                            |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+

{docstring_params_header_fine}

{docstring_params_table_header}
| unity_sum_coeff       | float         | 1e6                           | Coefficient of the term that enforces the unity sum constraint (:math:`\mathcal{{U}}`).                 |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| amp_penalty_max       | float         | 1e6                           | Maximum value of the regularization parameter that penalizes negative values of                       |
|                       |               |                               | the spectral function (:math:`\mathcal{{Q}}`).                                                          |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| amp_penalty_divisor   | float         | 100                           | Divisor used to reduce the regularization parameter that penalizes negative values of                 |
|                       |               |                               | the spectral function.                                                                                |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| der_penalty_init      | float         | 1e-3                          | Initial value of the regularization parameters that penalize large derivatives of the solution        |
|                       |               |                               | (:math:`\mathcal{{D}}`).                                                                                |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
| der_penalty_coeff     | float         | 2.0                           | Coefficient used to increase the regularization parameters that penalize large derivatives of         |
|                       |               |                               | the solution (:math:`f`).                                                                             |
+-----------------------+---------------+-------------------------------+-------------------------------------------------------------------------------------------------------+
""")

#
# SomCore.run()
#

c.add_method("void run(**som::accumulate_parameters_t)",
             doc = f"""
Accumulate particular solutions and compute the final solution using the standard SOM criterion
===============================================================================================

{docstring_accumulate}
""")

#
# SomCore.particular_solutions()
#

c.add_method(name = "particular_solutions",
             signature = "std::vector<std::pair<configuration, double>> get_particular_solutions (long i)",
             doc = """Accumulated particular solutions and their respective values of the objective function :math:`\chi^2` for the i-th diagonal matrix element.
The returned list includes only those solutions stored locally in the calling MPI process.""")

#
# SomCore.solution() and SomCore.solutions
#

c.add_method(name = "solution",
             signature = "configuration get_solution (long i)",
             doc = """Final solution for the i-th diagonal matrix element of the observable""")

c.add_property(name = "solutions",
               getter = cfunction("std::vector<configuration> get_solutions ()"),
               doc = """Final solutions, one per diagonal matrix element of the observable""")

#
# SomCore.objf() and SomCore.objf_list
#

c.add_method(name = "objf",
             signature = "double get_objf (long i)",
             doc = """Value of the objective function :math:`\chi^2` of the final solution for the i-th diagonal matrix element""")

c.add_property(name = "objf_list",
               getter = cfunction("std::vector<double> get_objf ()"),
               doc = """Values of the objective function :math:`\chi^2` of the final solutions, one value per a diagonal matrix element of the observable""")

#
# SomCore.histogram() and SomCore.histograms
#

c.add_method(name = "histogram",
             signature = "std::optional<histogram> get_histogram(long i)",
             doc = """Accumulated :math:`\chi` histogram for the i-th diagonal matrix element of the observable""")

c.add_property(name = "histograms",
               getter = cfunction(signature = "std::optional<std::vector<histogram>> get_histograms ()"),
               doc = """Accumulated :math:`\chi` histograms, one per diagonal matrix element of the observable""")

#
# Other attributes of SomCore
#

c.add_property(name = "observable_kind",
               getter = cfunction("observable_kind get_observable_kind()"),
               doc = """Kind of the observable""")

c.add_property(name = "dim",
               getter = cfunction("long get_dim()"),
               doc = """Matrix dimension of the observable""")

c.add_property(name = "last_accumulate_parameters",
               getter = cfunction("som::accumulate_parameters_t get_last_accumulate_parameters ()"),
               doc = """Set of parameters used in the last call to accumulate() """)

c.add_property(name = "accumulate_status",
               getter = cfunction("int get_accumulate_status ()"),
               doc = """Status of the accumulate() on exit """)

c.add_property(name = "objf_min",
               getter = cfunction("std::vector<double> get_objf_min ()"),
               doc = """Minimum of the objective function :math:`\chi^2` over all accumulated particular solutions (one value per diagonal matrix element of the observable)""")

module.add_class(c)

#
# fill_refreq()
#

module.add_function("void fill_refreq(gf_view<refreq> g_w, som_core cont, bool with_binning = false)",
                    doc = """Fill a real-frequency observable from a computed SOM solution""")

module.add_function("void fill_refreq(gf_view<refreq> g_w, observable_kind kind, std::vector<configuration> solutions, bool with_binning = false)",
                    doc = """Fill a real-frequency observable from a list of solutions (one solution per diagonal matrix element of the observable)""")

#
# compute_tail()
#

module.add_function("nda::array<dcomplex, 3> compute_tail(int max_order, som_core cont)",
                    doc = """Compute tail coefficients from a computed SOM solution""")

module.add_function("nda::array<dcomplex, 3> compute_tail(int max_order, observable_kind kind, std::vector<configuration> solutions)",
                    doc = """Compute tail coefficients from a list of solutions (one solution per diagonal matrix element of the observable)""")

#
# reconstruct()
#

module.add_function("void reconstruct(gf_view<imtime> g, som_core cont)",
                    doc = """Reconstruct an observable in the imaginary-time representation from a computed SOM solution""")
module.add_function("void reconstruct(gf_view<imfreq> g, som_core cont)",
                    doc = """Reconstruct an observable in the imaginary-frequency representation from a computed SOM solution""")
module.add_function("void reconstruct(gf_view<legendre> g, som_core cont)",
                    doc = """Reconstruct an observable in the Legendre polynomial basis from a computed SOM solution""")

module.add_function("void reconstruct(gf_view<imtime> g, observable_kind kind, std::vector<configuration> solutions)",
                    doc = """Reconstruct an observable in the imaginary-time representation from a list of solutions (one solution per diagonal matrix element of the observable)""")
module.add_function("void reconstruct(gf_view<imfreq> g, observable_kind kind, std::vector<configuration> solutions)",
                    doc = """Reconstruct an observable in the imaginary-frequency representation from a list of solutions (one solution per diagonal matrix element of the observable)""")
module.add_function("void reconstruct(gf_view<legendre> g, observable_kind kind, std::vector<configuration> solutions)",
                    doc = """Reconstruct an observable in the Legendre polynomial basis from a list of solutions (one solution per diagonal matrix element of the observable)""")

module.generate_code()
