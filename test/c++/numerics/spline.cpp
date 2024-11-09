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

#define NDA_ENFORCE_BOUNDCHECK
#include <som/numerics/spline.hpp>

using nda::vector;
using std::abs;
using std::copysign;
using std::exp;
using std::sin;

using namespace som;

// Function to interpolate
double f(double x) { return 2 * exp(-0.5 * (x - 1) * (x - 1)) * sin(0.5 * x); }

int const N = 20001;
double const xmin = -10, xmax = 10;
double const h = (xmax - xmin) / (N - 1);

// Regular mesh
double r_mesh(int i) { return xmin + i * h; }
// Irregular meshes
double ir_mesh1(int i) {
  auto x = r_mesh(i);
  return 0.6 * (x + copysign(3 - 3 * exp(-2 * abs(x)), x));
}
double ir_mesh2(int i) {
  auto x = r_mesh(i);
  return 0.6 * (x + copysign(1 - 1 * exp(-5 * abs(x)), x));
}

TEST(spline_test, spline) {
  vector<double> x(N), y(N);

  for(int i = 0; i < N; ++i) {
    x(i) = ir_mesh1(i);
    y(i) = f(x(i));
  }

  spline sp(x, y);

  vector<double> ref(N), result(N);
  for(int i = 0; i < N; ++i) {
    double z = r_mesh(i);
    ref(i) = f(z);
    result(i) = sp(z);
  }
  EXPECT_ARRAY_NEAR(result, ref, 1e-9);

  for(int i = 0; i < N; ++i) {
    double z = ir_mesh2(i);
    ref(i) = f(z);
    result(i) = sp(z);
  }
  EXPECT_ARRAY_NEAR(result, ref, 1e-9);
}

TEST(spline_test, regular_spline) {
  vector<double> y(N);

  for(int i = 0; i < N; ++i) y(i) = f(r_mesh(i));

  regular_spline sp(xmin, xmax, y);

  vector<double> ref(N), result(N);
  for(int i = 0; i < N; ++i) {
    double z = ir_mesh1(i);
    ref(i) = f(z);
    result(i) = sp(z);
  }
  EXPECT_ARRAY_NEAR(result, ref, 1e-13);

  for(int i = 0; i < N; ++i) {
    double z = ir_mesh2(i);
    ref(i) = f(z);
    result(i) = sp(z);
  }
  EXPECT_ARRAY_NEAR(result, ref, 1e-13);
}
