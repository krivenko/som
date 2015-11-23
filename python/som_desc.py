# Generated automatically using the command :
# c++2py.py ../c++/som.hpp -p -mpytriqs.applications.analytic_continuation.som -o som --moduledoc "The Stochastic Optimization Method"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.analytic_continuation.som", doc = "The Stochastic Optimization Method")

# All the triqs C++/Python modules

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/som.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
using namespace som;
""")

# The class continuation
c = class_(
        py_type = "Continuation",  # name of the python class
        c_type = "continuation",   # name of the C++ class
)

module.add_class(c)

module.generate_code()