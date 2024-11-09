/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko
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

#include <cmath>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <h5/h5.hpp>

#include <som/configuration.hpp>
#include <som/kernels/fermiongf_imtime.hpp>
#include <som/solution_worker.hpp>

// #define REPACKAGE_ARCHIVE

using namespace nda;
using namespace som;
using triqs::utility::clock_callback;

// Configurations are close
::testing::AssertionResult conf_are_close(configuration const &x,
                                          configuration const &y,
                                          double tolerance = 1.e-7) {
  if(x.size() != y.size())
    return ::testing::AssertionFailure()
           << "Comparing two configuration of different size "
           << "\n X = " << x << "\n Y = " << y;

  for(auto const &[rx, ry] : itertools::zip(x, y)) {
    if(std::abs(rx.center - ry.center) > tolerance)
      return ::testing::AssertionFailure() << "Rectangles " << rx << " and "
                                           << ry << " have different centers";
    if(std::abs(rx.width - ry.width) > tolerance)
      return ::testing::AssertionFailure() << "Rectangles " << rx << " and "
                                           << ry << " have different widths";
    if(std::abs(rx.height - ry.height) > tolerance)
      return ::testing::AssertionFailure() << "Rectangles " << rx << " and "
                                           << ry << " have different heights";
  }

  return ::testing::AssertionSuccess();
}

#define EXPECT_CONF_NEAR(X, ...) EXPECT_TRUE(conf_are_close(X, __VA_ARGS__))

struct solution_worker_test : public ::testing::Test {
protected:
  h5::file arch;

  double beta = {};
  triqs::mesh::imtime mesh;
  array<double, 1> g_tau;
  array<double, 1> error_bars_tau;
  matrix<double> cov_matrix_tau;
  double filtering_level = 0;

  using obj_function =
      objective_function<kernel<FermionGf, triqs::mesh::imtime>>;

  constexpr static char arch_open_mode =
#ifdef REPACKAGE_ARCHIVE
      'a';
#else
      'r';
#endif

public:
  solution_worker_test() : arch("solution_worker.ref.h5", arch_open_mode) {

    h5_read(arch, "beta", beta);
    h5_read(arch, "g_tau", g_tau);
    h5_read(arch, "error_bars_tau", error_bars_tau);
    h5_read(arch, "cov_matrix_tau", cov_matrix_tau);
    h5_read(arch, "filtering_level", filtering_level);

    mesh = {beta, triqs::mesh::Fermion, first_dim(g_tau)};
  }
};

TEST_F(solution_worker_test, RandomConfig) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, error_bars_tau);

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963160;
  params.t = 100;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  auto solution = worker(10);

#ifdef REPACKAGE_ARCHIVE
  h5_write(arch, "RandomConfig_output", solution);
#else
  configuration solution_ref(ci);
  h5_read(arch, "RandomConfig_output", solution_ref);

  EXPECT_CONF_NEAR(solution_ref, solution);
#endif
}

TEST_F(solution_worker_test, StartConfig) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, error_bars_tau);

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963160;
  params.t = 100;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  configuration init_config(ci);
  h5_read(arch, "StartConfig_input", init_config);

  auto solution = worker(init_config);

#ifdef REPACKAGE_ARCHIVE
  h5_write(arch, "StartConfig_output", solution);
#else
  configuration solution_ref(ci);
  h5_read(arch, "StartConfig_output", solution_ref);

  EXPECT_CONF_NEAR(solution_ref, solution);
#endif
}

TEST_F(solution_worker_test, RandomConfig_cov_matrix) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, cov_matrix_tau, filtering_level);

  nda::array<double, 1> sigma2_ref;
  h5_read(arch, "cov_matrix_tau_sigma2", sigma2_ref);
  EXPECT_ARRAY_NEAR(sigma2_ref, of.get_sigma2());

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963160;
  params.t = 100;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  auto solution = worker(10);

#ifdef REPACKAGE_ARCHIVE
  h5_write(arch, "RandomConfig_output_cov_matrix", solution);
#else
  configuration solution_ref(ci);
  h5_read(arch, "RandomConfig_output_cov_matrix", solution_ref);

  EXPECT_CONF_NEAR(solution_ref, solution);
#endif
}

TEST_F(solution_worker_test, StartConfig_cov_matrix) {
  cache_index ci;
  kernel<FermionGf, triqs::mesh::imtime> kern(mesh);
  obj_function of(kern, g_tau, cov_matrix_tau, filtering_level);

  nda::array<double, 1> sigma2_ref;
  h5_read(arch, "cov_matrix_tau_sigma2", sigma2_ref);
  EXPECT_ARRAY_NEAR(sigma2_ref, of.get_sigma2());

  auto params = worker_parameters_t({-3.0, 3.0});
  params.random_seed = 963160;
  params.t = 100;
  params.min_rect_width = 0.001;
  params.min_rect_weight = 0.001;

  solution_worker<kernel<FermionGf, triqs::mesh::imtime>> worker(
      of, 1.0, ci, params, clock_callback(-1), 10);

  configuration init_config(ci);
  h5_read(arch, "StartConfig_input", init_config);

  auto solution = worker(init_config);

#ifdef REPACKAGE_ARCHIVE
  h5_write(arch, "StartConfig_output_cov_matrix", solution);
#else
  configuration solution_ref(ci);
  h5_read(arch, "StartConfig_output_cov_matrix", solution_ref);

  EXPECT_CONF_NEAR(solution_ref, solution);
#endif
}
