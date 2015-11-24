#pragma once

#include <vector>

#include <triqs/arrays/vector.hpp>
#include <triqs/gfs.hpp>

#include "run_parameters.hpp"
#include "configuration.hpp"
#include "kernels.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

class som_core {

 // The right-hand side of the Fredholm integral equation
 // One vector per diagonal matrix element of the Green's function
 std::vector<vector<double>> rhs;
 // Error bars of the RHS, see eq. 30
 // One vector per diagonal matrix element of the Green's function
 std::vector<vector<double>> error_bars;

public:

 som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S);

 TRIQS_WRAP_ARG_AS_DICT // Wrap the parameters as a dictionary in python with the clang tool
 void run(run_parameters_t const& p);
};

}
