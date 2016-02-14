# Generated automatically using the command :
# c++2py.py ../c++/som_core.hpp -p -mpytriqs.applications.analytic_continuation.som -o som --moduledoc "The Stochastic Optimization Method"
from wrap_generator import *

# The module
module = module_(app_name="triqs_som", full_name = "som", doc = "The Stochastic Optimization Method")

# All the triqs C++/Python modules
module.use_module('gf', 'triqs')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/som_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
using namespace som;
#include "converters.hxx"
""")

# The class som_core
c = class_(
        py_type = "SomCore",  # name of the python class
        c_type = "som_core",   # name of the C++ class
)

c.add_constructor("""(gf_view<imtime> g_tau, gf_view<imtime> S_tau)""",
                  doc = """ """)

c.add_constructor("""(gf_view<imfreq> g_iw, gf_view<imfreq> S_iw)""",
                  doc = """ """)

c.add_constructor("""(gf_view<legendre> g_l, gf_view<legendre> S_l)""",
                  doc = """ """)

c.add_method("""void run (**som::run_parameters_t)""",
             doc = """  Parameter Name Type         Default Documentation

  max_rects      unsigned int 60      Maximum number of rectangles to represent spectra (K_{max}), should be below 70  """)

c.add_call(signature = "gf_view<refreq>(gf_view<refreq> g_w)")

module.add_class(c)

module.generate_code()