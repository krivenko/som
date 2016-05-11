from wrap_generator import *

# The module
module = module_(full_name = "histogram", doc = "Statistical histogram", app_name = "triqs_som")

# All the triqs C++/Python modules

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("<triqs/statistics/histograms.hpp>")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/arrays.hpp>
#include <triqs/python_tools/converters/h5.hpp>
using namespace triqs;
using namespace statistics;
""")

# The class histogram
c = class_(
        py_type = "Histogram",  # name of the python class
        c_type = "histogram",   # name of the C++ class
        c_type_absolute = "triqs::statistics::histogram",
        is_printable = True,
        hdf5 = True,
        arithmetic = "add_only"
)

c.add_constructor("""(int a, int b)""",
                  doc = """Constructor with mesh of integer values""")

c.add_constructor("""(double a, double b, long n_bins)""",
                  doc = """Constructor with mesh of double values""")

c.add_method("""double mesh_point (int n)""",
             doc = """Point on the mesh """)

c.add_property(name = "limits",
               getter = cfunction("std::pair<double,double> limits ()"),
               doc = """Returns boundaries of the histogram """)

c.add_property(name = "data",
               getter = cfunction("arrays::vector<double> data ()"),
               doc = """Read-only access to the data storage """)

c.add_property(name = "n_data_pts",
               getter = cfunction("uint64_t n_data_pts ()"),
               doc = """Norm of the stored data """)

c.add_property(name = "n_lost_pts",
               getter = cfunction("uint64_t n_lost_pts ()"),
               doc = """Number of discarded points """)

c.add_method("void clear ()", doc = """Reset all values to 0""")

c.add_len(doc = "Number of bins")

c.number_protocol['inplace_lshift'] = pyfunction(name ="__ilshift__", arity = 2)
c.number_protocol['inplace_lshift'].add_overload(calling_pattern = "<<", signature = "histogram(histogram h, double x)")

module.add_class(c)

module.add_function ("triqs::statistics::histogram pdf (triqs::statistics::histogram h)", doc = "Normalise histogram to get PDF")

module.add_function ("triqs::statistics::histogram cdf (triqs::statistics::histogram h)", doc = "Integrate and normalise histogram to get CDF")

module.generate_code()
