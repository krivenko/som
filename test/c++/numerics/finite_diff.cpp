/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#define NDA_ENFORCE_BOUNDCHECK
#include <som/numerics/finite_diff.hpp>

using namespace nda;
using namespace som;

// TODO: linspace()?
array<double, 1> const n = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
array<double, 1> const x = n + 0.1 * n * n;
array<double, 1> const f = exp(x / 5);

TEST(finite_diff_test, finite_diff_forward) {
  array<double, 1> result(f.size() - 1);
  finite_diff_forward(make_const_view(f), x, result());
  array<double, 1> ref = {
      0.2237061187158007,
      0.2846135935427021,
      0.37693190887020506,
      0.5196364339969416,
      0.7457025616026646,
      1.1139378090250882,
      1.7321497368002228, // NOLINT(modernize-use-std-numbers)
      2.8037481262723767,
      4.724126607940005,
      8.285770693825526};
  EXPECT_ARRAY_NEAR(result, ref, 1e-12);

  array<double, 1> dx(x.size() - 1);
  for(auto i : range(x.size() - 1)) dx(i) = x(i + 1) - x(i);

  finite_diff_forward_dx(make_const_view(f), dx, result());
  EXPECT_ARRAY_NEAR(result, ref, 1e-12);
}

TEST(finite_diff_test, finite_2_symm) {
  array<double, 1> result(f.size() - 2);
  finite_diff_2_symm(make_const_view(f), x, result());
  array<double, 1> ref = {0.05075622902241774,
                          0.06594165380535928,
                          0.08919032820421023,
                          0.12559229311429074,
                          0.18411762371121196,
                          0.2810054217159704,
                          0.44649932894673083,
                          0.7386071083337036,
                          1.2720157449591143};
  EXPECT_ARRAY_NEAR(result, ref, 1e-12);

  array<double, 1> dx(x.size() - 1);
  for(auto i : range(x.size() - 1)) dx(i) = x(i + 1) - x(i);

  finite_diff_2_symm_dx(make_const_view(f), dx, result());
  EXPECT_ARRAY_NEAR(result, ref, 1e-12);
}
