#include <triqs/test_tools/arrays.hpp>

#include "som_core.hpp"

using namespace triqs::gfs;
using namespace som;

const double beta = 10;

TEST(gf, gf_imtime) {

 auto g_tau = gf<imtime>{{beta, Fermion, 200}, {1,1}};
 auto s_tau = gf<imtime>{{beta, Fermion, 200}, {1,1}};

 auto g_w = gf<refreq>{{-5, 5, 1000}, {1,1}};

 som_core continuation(g_tau,s_tau);

 auto params = run_parameters_t{};
 continuation.run(params);

 //g_w << continuation;

}

MAKE_MAIN;
