/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2015 by I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <limits>

#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/gtsv.hpp>

namespace som {

using namespace triqs::arrays;

// Cubic spline interpolator
class spline {

 vector<double> x, y; // vectors of knot coordinates
 vector<double> a, b; // spline coefficients
 int s;               // size of x and y

public:

 spline(vector<double> const& x, vector<double> const& y) :
  x(x), y(y), s(x.size()), a(x.size()-1,0), b(x.size()-1,0) {
  if(s != y.size())
  TRIQS_RUNTIME_ERROR << "spline : inconsistent sizes of x and y vectors in construction, " << s << " vs " << y.size();
  if(s < 2) TRIQS_RUNTIME_ERROR << "spline : need at least 2 knots";

  double max_x = -std::numeric_limits<double>::infinity();
  for(auto x_ : x) {
   if(x_ < max_x) TRIQS_RUNTIME_ERROR << "spline : vector x is not sorted in the ascending order";
   max_x = x_;
  }

  vector<double> ld(s-1), d(s), ud(s-1), k(s);
  d(0) = 2/(x(1) - x(0));
  k(0) = 3*(y(1) - y(0))/((x(1)-x(0))*(x(1)-x(0)));
  for(int i = 1; i < s-1; ++i) {
   double xm = x(i-1), xi = x(i), xp = x(i+1);
   d(i) = 2*(1.0/(xi - xm) + 1.0/(xp - xi));
   ld(i-1) = ud(i-1) = 1/(xi - xm);
   k(i) = 3*((y(i+1) - y(i))/((xp-xi)*(xp-xi)) + (y(i) - y(i-1))/((xi-xm)*(xi-xm)));
  }
  d(s-1) = 2/(x(s-1) - x(s-2));
  ld(s-2) = ud(s-2) = 1/(x(s-1) - x(s-2));
  k(s-1) = 3*(y(s-1) - y(s-2))/((x(s-1)-x(s-2))*(x(s-1)-x(s-2)));

  blas::gtsv(ld,d,ud,k);
  for(int i = 0; i < s-1; ++i) {
   a(i) = k(i)*(x(i+1) - x(i)) - (y(i+1) - y(i));
   b(i) = -k(i+1)*(x(i+1) - x(i)) + (y(i+1) - y(i));
  }
 }

 inline double operator()(double z) const {
  if(z <= x(0)) return y(0);
  else if(z >= x(s-1)) return y(s-1);
  else {
   int left = 0, right = s-1;
   while(left != right-1) {
    int mid = (left+right)/2;
    double mid_x = x(mid);
    if(z < mid_x) right = mid;
    else if(z > mid_x) left = mid;
    else break;
   }
   double t = (z - x(left))/(x(right) - x(left));
   double tt = 1 - t;
   return tt*y(left)  + t*y(right) + t*tt*(a(left)*tt + b(left)*t);
  }
 }
};

// Cubic spline interpolator on a regular mesh
class regular_spline {
 // TODO
};

}