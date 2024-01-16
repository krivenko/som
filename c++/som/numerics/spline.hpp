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
#pragma once

#include <nda/nda.hpp>

namespace som {

// Cubic spline interpolator
class spline {

  nda::array<double, 1> x, y; // vectors of knot coordinates
  nda::array<double, 1> a, b; // spline coefficients
  size_t s = 0;               // size of x and y

public:
  spline() = default;
  spline(nda::vector<double> const& x, nda::vector<double> const& y);

  double operator()(double z) const;
};

// Cubic spline interpolator on a regular mesh
class regular_spline {

  double x_min = 0, x_max = 0; // abscissae of the leftmost and rightmost knots
  nda::array<double, 1> y;     // vector of knot ordinates
  nda::array<double, 1> a, b;  // spline coefficients
  size_t s = 0;                // size of y
  double dx = 0;               // mesh step

public:
  regular_spline() = default;
  regular_spline(double x_min, double x_max, nda::vector<double> const& y);

  double operator()(double z) const;
};

} // namespace som
