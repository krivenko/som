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
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>

#include <cmath>

#include <som/numerics/simpson.hpp>

using namespace som;
using namespace std;

TEST(simpson_test, exp_sin) {
  auto f = [](double x) { return exp(2 * x) * sin(3 * x); };

  double ref = exp(6) * (2 * sin(9) - 3 * cos(9)) / 13;
  ref -= exp(2) * (2 * sin(3) - 3 * cos(3)) / 13;
  EXPECT_NEAR(ref, adaptive_simpson(f, 1.0, 3.0, 1e-10), 1e-10);
}

TEST(simpson_test, sin_inv_x) {
  auto f = [](double x) { return sin(1 / x) / (1 + 2 * x + x * x); };
  EXPECT_NEAR(0.162985567, adaptive_simpson(f, 0.01, 1.0, 1e-10), 1e-9);
}

TEST(simpson_test, primitive) {
  auto f = [](double x) { return x * x; };
  nda::array<double, 1> ref(11);

  nda::for_each(ref.shape(), [&ref](int i) {
    double y = 0.1 + 0.1 * i;
    ref(i) = y * y * y / 3 - 0.001 / 3;
  });
  EXPECT_ARRAY_NEAR(ref, primitive(f, 0.1, 1.1, 11, 1e-10, false), 1e-10);

  nda::for_each(ref.shape(), [&ref](int i) {
    double y = -1.1 + 0.1 * i;
    ref(i) = y * y * y / 3 + 0.001 / 3;
  });
  EXPECT_ARRAY_NEAR(ref, primitive(f, -1.1, -0.1, 11, 1e-10, true), 1e-10);
}
