#pragma once

#include <vector>

#include <triqs/arrays/vector.hpp>
#include <triqs/gfs.hpp>
#include <triqs/utility/variant.hpp>

#include "run_parameters.hpp"
#include "configuration.hpp"
#include "kernels.hpp"
#include "objective_function.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;
using triqs::utility::variant;

class som_core {

 using vector_of_real_or_complex_vector_t =
 variant<std::vector<vector<double>>,std::vector<vector<double>>>;

 // The right-hand side of the Fredholm integral equation
 // One vector per diagonal matrix element of the Green's function
 vector_of_real_or_complex_vector_t rhs;
 // Error bars of the RHS, see eq. 30
 // One vector per diagonal matrix element of the Green's function
 vector_of_real_or_complex_vector_t error_bars;

 // Resulting configurations
 std::vector<configuration> results;

public:

 som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S);
 //som_core(gf_const_view<imfreq> g_iw, gf_const_view<imfreq> S);

 TRIQS_WRAP_ARG_AS_DICT // Wrap the parameters as a dictionary in python with the clang tool
 void run(run_parameters_t const& p);

 // Fill gf<refreq> with obtained results
 friend void triqs_gf_view_assign_delegation(gf_view<refreq> g_w, som_core const& cont) {
  // TODO
 }

};

}
