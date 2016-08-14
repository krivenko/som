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
#include <triqs/h5.hpp>

#include "kernels/fermiongf_imtime.hpp"
#include "solution_worker.hpp"

using namespace som;
using namespace triqs::arrays;
using triqs::utility::clock_callback;

double beta;
gf_mesh<imtime> mesh;
array<double,1> g_tau, s_tau;

struct Env : public ::testing::Environment {
 virtual void SetUp() {
  triqs::h5::file arch("solution_worker.ref.h5", H5F_ACC_RDONLY);

  h5_read(arch, "beta", beta);
  h5_read(arch, "g_tau", g_tau);
  h5_read(arch, "s_tau", s_tau);

  mesh = {beta, Fermion, static_cast<int>(first_dim(g_tau))};
 }
};
::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(new Env);

using obj_function = objective_function<kernel<FermionGf,imtime>>;

TEST(solution_worker,RandomConfig) {
 cache_index ci;
 kernel<FermionGf,imtime> kern(mesh);
 obj_function of(kern, g_tau, s_tau);

 auto params = run_parameters_t({-3.0,3.0});
 params.random_seed = 963162;
 params.t = 100;
 params.min_rect_width = 0.001;
 params.min_rect_weight = 0.001;

 solution_worker<kernel<FermionGf,imtime>> worker(of, 1.0, ci, params, clock_callback(-1), 10);

 auto solution = worker(10);

 configuration solution_ref(ci);
 triqs::h5::file arch("solution_worker.ref.h5", H5F_ACC_RDONLY);
 h5_read(arch, "RandomConfig_output", solution_ref, ci);

 EXPECT_EQ(solution_ref, solution);
}

TEST(solution_worker,StartConfig) {
 cache_index ci;
 kernel<FermionGf,imtime> kern(mesh);
 obj_function of(kern, g_tau, s_tau);

 auto params = run_parameters_t({-3.0,3.0});
 params.random_seed = 963162;
 params.t = 100;
 params.min_rect_width = 0.001;
 params.min_rect_weight = 0.001;

 solution_worker<kernel<FermionGf,imtime>> worker(of, 1.0, ci, params, clock_callback(-1), 10);

 triqs::h5::file arch("solution_worker.ref.h5", H5F_ACC_RDONLY);
 configuration init_config(ci);
 h5_read(arch, "StartConfig_input", init_config, ci);

 auto solution = worker(init_config);

 configuration solution_ref(ci);
 h5_read(arch, "StartConfig_output", solution_ref, ci);

 EXPECT_EQ(solution_ref, solution);
}

MAKE_MAIN;
