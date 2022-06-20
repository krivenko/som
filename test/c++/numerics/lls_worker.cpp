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

#include <cmath>
#include <stdexcept>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <triqs/utility/exceptions.hpp>

#include <som/numerics/lls_worker.hpp>

using namespace som;

TEST(lls_worker, M_exceeds_N) {
  // clang-format off
  nda::matrix<double, nda::F_layout> A = {{1, 1, 2, 3},
                                          {1, 2, 3, 4},
                                          {1, 3, 4, 5},
                                          {1, 4, 5, 6},
                                          {1, 5, 6, 7}};
  // clang-format on
  auto M = int(first_dim(A));
  auto N = int(second_dim(A));

  nda::vector<double> b = {7, 7, 8, 8, 9};
  nda::vector<double> x(N);

  auto worker = lls_worker(M, N);
  auto rank = worker(A, b, x);

  EXPECT_EQ(rank, 2);
  nda::vector<double> ref = {29.0 / 15, -53.0 / 30, 1.0 / 6, 2.1};
  EXPECT_ARRAY_NEAR(x, ref, 1e-12);
  nda::array<double, 1> sigma_ref = {std::sqrt(7.5 * (19 + std::sqrt(353))),
                                     std::sqrt(7.5 * (19 - std::sqrt(353))),
                                     0,
                                     0};
  EXPECT_ARRAY_NEAR(worker.singular_values(), sigma_ref, 1e-12);
}

TEST(lls_worker, N_exceeds_M) {
  // clang-format off
  nda::matrix<double, nda::F_layout> A = {{1, 1, 1, 1, 1},
                                          {1, 2, 3, 4, 5},
                                          {2, 3, 4, 5, 6},
                                          {3, 4, 5, 6, 7}};
  // clang-format on
  auto M = int(first_dim(A));
  auto N = int(second_dim(A));

  nda::vector<double> b = {7, 7, 8, 8};
  nda::vector<double> x(N);

  auto worker = lls_worker(M, N);
  auto rank = worker(A, b, x);

  EXPECT_EQ(rank, 2);
  nda::vector<double> ref = {17.0 / 15, 5.0 / 6, 8.0 / 15, 7.0 / 30, -1.0 / 15};
  EXPECT_ARRAY_NEAR(x, ref, 1e-12);
  nda::array<double, 1> sigma_ref = {std::sqrt(7.5 * (19 + std::sqrt(353))),
                                     std::sqrt(7.5 * (19 - std::sqrt(353))),
                                     0,
                                     0};
  EXPECT_ARRAY_NEAR(worker.singular_values(), sigma_ref, 1e-12);
}

TEST(lls_worker, exceptions) {
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
  EXPECT_THROW(lls_worker(0, 10), triqs::runtime_error);
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-goto)
  EXPECT_THROW(lls_worker(10, 0), triqs::runtime_error);
}
