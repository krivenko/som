/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include <h5/h5.hpp>
#include <nda/gtest_tools.hpp>

#include <som/kernels/fermiongf_imtime.hpp>
#include <som/solution_worker.hpp>

using namespace nda;
using namespace som;
using triqs::utility::clock_callback;

struct solution_worker_test : public ::testing::Test {
protected:
  h5::file arch;

  double beta = {};
  triqs::mesh::imtime mesh;
  array<double, 1> g_tau;
  array<double, 1> error_bars_tau;
  matrix<double> cov_matrix_tau;
  double filtration_level = 0;

  using obj_function =
      objective_function<kernel<FermionGf, triqs::mesh::imtime>>;

public:
  solution_worker_test() : arch("solution_worker.ref.h5", 'r') {

    h5_read(arch, "beta", beta);
    h5_read(arch, "g_tau", g_tau);
    h5_read(arch, "error_bars_tau", error_bars_tau);
    h5_read(arch, "cov_matrix_tau", cov_matrix_tau);
    h5_read(arch, "filtration_level", filtration_level);

    mesh = {beta, triqs::mesh::Fermion, first_dim(g_tau)};
  }
};

TEST_F(solution_worker_test, RandomConfig) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, error_bars_tau);

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963162;
  params.t = 100;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  auto solution = worker(10);

  configuration solution_ref(ci);
  h5_read(arch, "RandomConfig_output", solution_ref);

  EXPECT_EQ(solution_ref, solution);
}

TEST_F(solution_worker_test, StartConfig) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, error_bars_tau);

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963162;
  params.t = 100;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  configuration init_config(ci);
  h5_read(arch, "StartConfig_input", init_config);

  auto solution = worker(init_config);

  configuration solution_ref(ci);
  h5_read(arch, "StartConfig_output", solution_ref);

  EXPECT_EQ(solution_ref, solution);
}

TEST_F(solution_worker_test, RandomConfig_CC) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, error_bars_tau);

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963162;
  params.t = 1000;
  params.t1 = 800;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;
  params.cc_update = true;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  auto solution = worker(10);

  configuration solution_ref(ci);
  h5_read(arch, "RandomConfig_output_CC", solution_ref);

  EXPECT_EQ(solution_ref, solution);
}

TEST_F(solution_worker_test, StartConfig_CC) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, error_bars_tau);

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963162;
  params.t = 1000;
  params.t1 = 800;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;
  params.cc_update = true;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  configuration init_config(ci);
  h5_read(arch, "StartConfig_input", init_config);

  auto solution = worker(init_config);

  configuration solution_ref(ci);
  h5_read(arch, "StartConfig_output_CC", solution_ref);

  EXPECT_EQ(solution_ref, solution);
}

TEST_F(solution_worker_test, RandomConfig_CC_cov_matrix) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, cov_matrix_tau, filtration_level);

  nda::array<double, 1> sigma2_ref;
  h5_read(arch, "cov_matrix_tau_sigma2", sigma2_ref);
  EXPECT_ARRAY_NEAR(sigma2_ref, of.get_sigma2());

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963162;
  params.t = 1000;
  params.t1 = 800;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;
  params.cc_update = true;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  auto solution = worker(10);

  configuration solution_ref(ci);
  h5_read(arch, "RandomConfig_output_CC_cov_matrix", solution_ref);

  EXPECT_EQ(solution_ref, solution);
}

TEST_F(solution_worker_test, StartConfig_CC_cov_matrix) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, cov_matrix_tau, filtration_level);

  nda::array<double, 1> sigma2_ref;
  h5_read(arch, "cov_matrix_tau_sigma2", sigma2_ref);
  EXPECT_ARRAY_NEAR(sigma2_ref, of.get_sigma2());

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963162;
  params.t = 1000;
  params.t1 = 800;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;
  params.cc_update = true;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  configuration init_config(ci);
  h5_read(arch, "StartConfig_input", init_config);

  auto solution = worker(init_config);

  configuration solution_ref(ci);
  h5_read(arch, "StartConfig_output_CC_cov_matrix", solution_ref);

  EXPECT_EQ(solution_ref, solution);
}
