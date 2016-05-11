# Generated automatically using the command :
# c++2py.py ../c++/som_core.hpp --compiler_options=-DCACHE_SIZE=0x10000 -p -m som -o som --appname triqs_som --moduledoc "The Stochastic Optimization Method"
from wrap_generator import *

# The module
module = module_(full_name = "core", doc = "The Stochastic Optimization Method", app_name = "triqs_som")

# All the triqs C++/Python modules
module.use_module('gf', 'triqs')
module.use_module('histogram', 'triqs_som')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("som_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/vector.hpp>
using namespace triqs::gfs;
using namespace som;
#include "./converters.hxx"
""")

module.add_enum("observable_kind", ["FermionGf","Susceptibility","Conductivity"], "som", "Kinds of observables")

# The class som_core
c = class_(
        py_type = "SomCore",  # name of the python class
        c_type = "som_core",   # name of the C++ class
        doc = r"",   # doc of the C++ class
)

c.add_constructor("""(gf_view<imtime> g_tau, gf_view<imtime> S_tau, som::observable_kind kind = FermionGf, double norm = 1.0)""",
                  doc = """Construct on imaginary-time quantities """)

c.add_constructor("""(gf_view<imfreq> g_iw, gf_view<imfreq> S_iw, som::observable_kind kind = FermionGf, double norm = 1.0)""",
                  doc = """Construct on imaginary-frequency quantities """)

c.add_constructor("""(gf_view<legendre> g_l, gf_view<legendre> S_l, som::observable_kind kind = FermionGf, double norm = 1.0)""",
                  doc = """Construct on quantities in Legendre polynomial basis """)

c.add_method("""void run (**som::run_parameters_t)""",
             doc = """+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| Parameter Name         | Type            | Default                       | Documentation                                                                                                                                 |
+========================+=================+===============================+===============================================================================================================================================+
| energy_window          | (double,double) |                               | Estimated lower and upper bounds of the spectrum Negative values of the lower bound will be reset to 0 for susceptibilities and conductivity  |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| random_seed            | int             | 34788 + 928374 * MPI.rank     | Seed for random number generator                                                                                                              |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| random_name            | str             | ""                            | Name of random number generator                                                                                                               |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| max_rects              | unsigned int    | 60                            | Maximum number of rectangles to represent spectra (:math:`K_{max}`), should be below 70                                                       |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| min_rect_width         | double          | 1e-3                          | Minimal width of a rectangle, in units of the energy window width                                                                             |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| min_rect_weight        | double          | 1e-3                          | Minimal weight of a rectangle, in units of the requested norm of a solution                                                                   |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| n_elementary_updates   | int             | 1000                          | Number of elementary updates per global update (:math:`T`)                                                                                    |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| distrib_d_max          | double          | 2                             | Maximal parameter of the power-law distribution function for the Metropolis algorithm                                                         |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| gamma                  | double          | 2                             | Proposal probability parameter :math:`\gamma`                                                                                                 |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| n_global_updates       | unsigned int    | 100                           | Number of global updates (:math:`F`)                                                                                                          |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| adjust_ngu             | bool            | true                          | Adjust the number of global updates automatically If `true`, use n_global_updates as a starting value                                         |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| adjust_ngu_n_solutions | unsigned int    | 10                            | Number of particular solutions used to adjust :math:`F`                                                                                       |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| adjust_ngu_kappa       | double          | 0.25                          | Limiting value of :math:`\kappa` used to adjust :math:`F`                                                                                     |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| n_solutions            | unsigned int    | 100                           | Number of particular solutions used in the final accumulation (:math:`L`)                                                                     |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| adjust_nsol            | bool            | true                          | Adjust the number of solutions used in the final accumulation If `true`, use n_solutions as a starting value                                  |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| adjust_nsol_good_d     | double          | 2.0                           | Maximal ratio :math:`D/D_\mathrm{min}` for a particular solution to be considered good                                                        |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| adjust_nsol_verygood_d | double          | 4.0/3.0                       | Maximal ratio :math:`D/D_\mathrm{min}` for a particular solution to be considered very good                                                   |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| adjust_nsol_ratio      | double          | 0.95                          | Critical ratio :math:`N_\mathrm{very good}/N_\mathrm{good}` to stop :math:`L`-adjustment procedure.                                           |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| make_histograms        | bool            | false                         | Accumulate histograms of objective function values                                                                                            |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| hist_max               | double          | 2.0                           | Right boundary of the histograms, in units of :math:`D_\mathrm{min}` (left boundary is always set to :math:`D_\mathrm{min}`)                  |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| hist_n_bins            | bool            | 100                           | Number of bins for the histograms                                                                                                             |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| verbosity              | int             | 3 on MPI rank 0, 0 otherwise. | Verbosity level                                                                                                                               |
+------------------------+-----------------+-------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+ """)

c.add_call(signature = "void(gf_view<refreq> g_w)", calling_pattern = "triqs_gf_view_assign_delegation(g_w, self_c)")
c.add_call(signature = "void(gf_view<imtime> g_tau)", calling_pattern = "triqs_gf_view_assign_delegation(g_tau, self_c)")
c.add_call(signature = "void(gf_view<imfreq> g_iw)", calling_pattern = "triqs_gf_view_assign_delegation(g_iw, self_c)")
# FIXME
#c.add_call(signature = "void(gf_view<legendre> g_l)", calling_pattern = "triqs_gf_view_assign_delegation(g_l, self_c)")

c.add_property(name = "last_run_parameters",
               getter = cfunction("som::run_parameters_t get_last_run_parameters ()"),
               doc = """Set of parameters used in the last call to run() """)

c.add_property(name = "histograms",
               getter = cfunction("std::vector<histogram> get_histograms ()"),
               doc = """Accumulated objective function histograms""")

module.add_class(c)

module.generate_code()
