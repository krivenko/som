/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016 by I. Krivenko
 *
 * SOM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * SOM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * SOM. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
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
