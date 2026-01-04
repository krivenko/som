/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2026 Igor Krivenko
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

#include <complex>
#include <stdexcept>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <triqs/utility/exceptions.hpp>

#include <som/numerics/ecqp_worker.hpp>

using namespace som;

TEST(ecqp_worker, solution) {
  // clang-format off
  nda::matrix<double, nda::F_layout> Q = {{1, 1, 1, 1},
                                          {1, 2, 1, 1},
                                          {1, 1, 3, 1},
                                          {1, 1, 1, 4}};
  nda::vector<double> f = {1, 1, 2, 2};
  nda::vector<double> x(4);
  // clang-format on

  // No constraints
  {
    nda::matrix<double, nda::F_layout> L(0, 4);
    nda::vector<double> d = {};

    auto worker = ecqp_worker(4, 0);
    double residue = worker(Q, f, L, d, x);

    nda::vector<double> x_ref = {0.166666666666, 0, 0.5, 0.333333333333};
    EXPECT_ARRAY_NEAR(x, x_ref, 1e-12);
    EXPECT_NEAR(residue, -1.833333333333, 1e-12);
  }
  // One constraint
  {
    nda::matrix<double, nda::F_layout> L = {{1, 1, 1, 1}};
    nda::vector<double> d = {0.3};

    auto worker = ecqp_worker(4, 1);
    double residue = worker(Q, f, L, d, x);

    nda::vector<double> x_ref = {-0.533333333333, 0, 0.5, 0.333333333333};
    EXPECT_ARRAY_NEAR(x, x_ref, 1e-12);
    EXPECT_NEAR(residue, -1.343333333333, 1e-12);
  }
  // Two constraints
  {
    nda::matrix<double, nda::F_layout> L = {{1, 1, 1, 1}, {1, 0, 0, 1}};
    nda::vector<double> d = {0.3, 0};

    auto worker = ecqp_worker(4, 2);
    double residue = worker(Q, f, L, d, x);

    nda::vector<double> x_ref = {-0.333333333333,
                                 -0.133333333333,
                                 // NOLINTNEXTLINE(modernize-use-std-numbers)
                                 0.433333333333,
                                 0.333333333333};
    EXPECT_ARRAY_NEAR(x, x_ref, 1e-12);
    EXPECT_NEAR(residue, -1.316666666666, 1e-12);
  }
}

TEST(ecqp_worker, exceptions) {
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
  EXPECT_THROW(ecqp_worker(0, 10), triqs::runtime_error);
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
  EXPECT_THROW(ecqp_worker(0, 10), triqs::runtime_error);

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
  EXPECT_THROW(ecqp_worker(1, -1), triqs::runtime_error);
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
  EXPECT_THROW(ecqp_worker(1, -1), triqs::runtime_error);

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
  EXPECT_THROW(ecqp_worker(5, 10), triqs::runtime_error);
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
  EXPECT_THROW(ecqp_worker(5, 10), triqs::runtime_error);

  nda::vector<double> f = {1, 1, 2, 2};
  nda::matrix<double, nda::F_layout> L0(0, 4);
  nda::vector<double> d0 = {};
  nda::vector<double> x(4);

  // Indefinite matrix
  {
    // clang-format off
    nda::matrix<double, nda::F_layout> Q = {{1, 1, 1, -1},
                                            {1, 2, 1, 1},
                                            {1, 1, 3, 1},
                                            {-1, 1, 1, 4}};
    // clang-format on

    auto worker = ecqp_worker(4, 0);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
    EXPECT_THROW(worker(Q, f, L0, d0, x), triqs::runtime_error);
  }

  // Rank-deficient L
  {
    // clang-format off
    nda::matrix<double, nda::F_layout> Q = {{1, 1, 1, 1},
                                            {1, 2, 1, 1},
                                            {1, 1, 3, 1},
                                            {1, 1, 1, 4}};
    // clang-format on
    nda::matrix<double, nda::F_layout> L = {{1, 1, 1, 1}, {1, 1, 1, 1}};
    nda::vector<double> d1 = {2, 2};

    auto worker = ecqp_worker(4, 2);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
    EXPECT_THROW(worker(Q, f, L, d1, x), triqs::runtime_error);
  }
}
