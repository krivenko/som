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

#include <cassert>
#include <cmath>
#include <limits>

#include <nda/linalg.hpp>
#include <triqs/utility/exceptions.hpp>

#include "spline.hpp"

namespace som {

using namespace nda;

////////////
// spline //
////////////

spline::spline(vector<double> const& x, vector<double> const& y)
   : x(x), y(y), a(x.size() - 1, 0), b(x.size() - 1, 0), s(x.size()) {
  if(s != y.size())
    TRIQS_RUNTIME_ERROR
        << "spline : inconsistent sizes of x and y vectors in construction, "
        << s << " vs " << y.size();
  if(s < 2) TRIQS_RUNTIME_ERROR << "spline : need at least 2 knots";

  double x_max = -std::numeric_limits<double>::infinity();
  for(auto x_ : x) {
    if(x_ < x_max)
      TRIQS_RUNTIME_ERROR
          << "spline : vector x is not sorted in the ascending order";
    x_max = x_;
  }

  vector<double> ld(s - 1), d(s), ud(s - 1), k(s);
  d(0) = 2 / (x(1) - x(0));
  k(0) = 3 * (y(1) - y(0)) / ((x(1) - x(0)) * (x(1) - x(0)));
  for(size_t i = 1; i < s - 1; ++i) {
    double xm = x(i - 1), xi = x(i), xp = x(i + 1);
    d(i) = 2 * (1.0 / (xi - xm) + 1.0 / (xp - xi));
    ld(i - 1) = ud(i - 1) = 1 / (xi - xm);
    k(i) = 3 * ((y(i + 1) - y(i)) / ((xp - xi) * (xp - xi)) +
                (y(i) - y(i - 1)) / ((xi - xm) * (xi - xm)));
  }
  d(s - 1) = 2 / (x(s - 1) - x(s - 2));
  ld(s - 2) = ud(s - 2) = 1 / (x(s - 1) - x(s - 2));
  k(s - 1) = 3 * (y(s - 1) - y(s - 2)) /
             ((x(s - 1) - x(s - 2)) * (x(s - 1) - x(s - 2)));

  if(lapack::gtsv(ld, d, ud, k))
    TRIQS_RUNTIME_ERROR << "spline: lapack::gtsv() has failed";
  for(size_t i = 0; i < s - 1; ++i) {
    a(i) = k(i) * (x(i + 1) - x(i)) - (y(i + 1) - y(i));
    b(i) = -k(i + 1) * (x(i + 1) - x(i)) + (y(i + 1) - y(i));
  }
}

double spline::operator()(double z) const {
  assert(std::isfinite(z));
  if(z <= x(0))
    return y(0);
  else if(z >= x(s - 1))
    return y(s - 1);
  else {
    size_t left = 0, right = s - 1;
    while(left != right - 1) {
      size_t mid = (left + right) / 2;
      double mid_x = x(mid);
      if(z < mid_x)
        right = mid;
      else if(z > mid_x)
        left = mid;
      else
        break;
    }
    double t = (z - x(left)) / (x(right) - x(left));
    double tt = 1 - t;
    return tt * y(left) + t * y(right) + t * tt * (a(left) * tt + b(left) * t);
  }
}

////////////////////
// regular_spline //
////////////////////

regular_spline::regular_spline(double x_min,
                               double x_max,
                               vector<double> const& y)
   : x_min(x_min)
   , x_max(x_max)
   , y(y)
   , a(y.size() - 1, 0)
   , b(y.size() - 1, 0)
   , s(y.size())
   , dx((x_max - x_min) / double(y.size() - 1)) {
  if(s < 2) TRIQS_RUNTIME_ERROR << "regular_spline : need at least 2 knots";
  if(x_min == x_max) {
    TRIQS_RUNTIME_ERROR
        << "regular_spline : x_min and x_max must not coincide ";
  }

  array<double, 1> ld = (1 / dx) * ones<double>(s - 1);
  array<double, 1> d = (4 / dx) * ones<double>(s);
  array<double, 1> ud = (1 / dx) * ones<double>(s - 1);

  array<double, 1> k(s);
  d(0) = d(s - 1) = 2 / dx;

  k(0) = 3 * (y(1) - y(0)) / (dx * dx);
  for(int i = 1; i < s - 1; ++i) k(i) = 3 * (y(i + 1) - y(i - 1)) / (dx * dx);
  k(s - 1) = 3 * (y(s - 1) - y(s - 2)) / (dx * dx);

  if(lapack::gtsv(ld, d, ud, k))
    TRIQS_RUNTIME_ERROR << "spline: lapack::gtsv() has failed";
  for(int i = 0; i < s - 1; ++i) {
    a(i) = k(i) * dx - (y(i + 1) - y(i));
    b(i) = -k(i + 1) * dx + (y(i + 1) - y(i));
  }
}

double regular_spline::operator()(double z) const {
  assert(std::isfinite(z));
  if(z <= x_min)
    return y(0);
  else if(z >= x_max)
    return y(s - 1);
  else {
    int left = std::floor((z - x_min) / dx);
    int right = left + 1;
    double t = (z - (x_min + left * dx)) / dx;
    double tt = 1 - t;
    return tt * y(left) + t * y(right) + t * tt * (a(left) * tt + b(left) * t);
  }
}

} // namespace som
