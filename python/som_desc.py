# Generated automatically using the command :
# c++2py.py ../c++/som_core.hpp -p -mpytriqs.applications.analytic_continuation.som -o som --moduledoc "The Stochastic Optimization Method"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.analytic_continuation.som", doc = "The Stochastic Optimization Method")

# All the triqs C++/Python modules
module.use_module('gf')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/som_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
using namespace triqs::gfs;
using namespace som;
#include "./converters.hxx"
""")

# The class som_core
c = class_(
        py_type = "SomCore",  # name of the python class
        c_type = "som_core",   # name of the C++ class
)

c.add_constructor("""(gf_view<imtime> g_tau, gf_view<imtime> S)""",
                  doc = """ """)

c.add_method("""void run (**som::run_parameters_t)""",
             doc = """  Parameter Name Type         Default Documentation

  max_rects      unsigned int 60      Maximum number of rectangles to represent spectra (K_{max}) Should be below 70  """)

module.add_class(c)

module.generate_code()