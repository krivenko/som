#pragma once

#include <triqs/gfs.hpp>
#include "run_parameters.hpp"

namespace som {

using namespace triqs::gfs;

class som_core {

public:

 som_core(gf_const_view<imtime> g_tau);

 TRIQS_WRAP_ARG_AS_DICT // Wrap the parameters as a dictionary in python with the clang tool
 void run(run_parameters_t const& p);
};

}