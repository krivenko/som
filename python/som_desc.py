# Generated automatically using the command :
# c++2py.py ../c++/som_core.hpp --compiler_options=-DCACHE_SIZE=0xffff -p -m som -o som --appname triqs_som --moduledoc "The Stochastic Optimization Method"
from wrap_generator import *

# The module
module = module_(full_name = "som", doc = "The Stochastic Optimization Method", app_name = "triqs_som")

# All the triqs C++/Python modules
module.use_module('gf', 'triqs')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/som_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
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
             doc = """+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| Parameter Name       | Type                      | Default                   | Documentation                                                                                                                                 |
+======================+===========================+===========================+===============================================================================================================================================+
| energy_window        | std::pair<double, double> |                           | Estimated lower and upper bounds of the spectrum Negative values of the lower bound will be reset to 0 for susceptibilities and conductivity  |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| random_seed          | int                       | 34788 + 928374 * MPI.rank | Seed for random number generator                                                                                                              |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| random_name          | str                       | ""                        | Name of random number generator                                                                                                               |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| max_rects            | unsigned int              | 60                        | Maximum number of rectangles to represent spectra (K_{max}), should be below 70                                                               |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| min_rect_width       | double                    | 1e-3                      | Minimal width of a rectangle, in units of the energy window width                                                                             |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| min_rect_weight      | double                    | 1e-3                      | Minimal weight of a rectangle, in units of the requested norm of a solution                                                                   |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| n_elementary_updates | int                       | 1000                      | Number of elementary updates per global update (T)                                                                                            |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| distrib_d_max        | double                    | 2                         | Maximal parameter of the power-law distribution function for the Metropolis algorithm                                                         |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| gamma                | double                    | 2                         | Proposal probability parameter :math:`\gamma`                                                                                                 |
+----------------------+---------------------------+---------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+ """)

c.add_call(signature = "gf_view<refreq>(gf_view<refreq> g_w)")

module.add_class(c)

module.generate_code()
