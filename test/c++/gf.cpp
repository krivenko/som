#include <triqs/test_tools/arrays.hpp>

#include "som_core.hpp"

using namespace triqs::gfs;
using namespace som;

const double beta = 10;

TEST(gf, gf_imtime) {

 auto g_tau = gf<imtime>{{beta, Fermion, 200}, {2,2}};
 auto s_tau = gf<imtime>{{beta, Fermion, 200}, {2,2}};

 auto g_w = gf<refreq>{{-5, 5, 1000}, {2,2}};

 som_core continuation(g_tau,s_tau);

 auto params = run_parameters_t({-5,5});
 //continuation.run(params);

 g_w() = continuation;

 triqs::h5::file G_file("gf.out.h5",'w');
 h5_write(G_file,"g_w",g_w);
}

MAKE_MAIN;
