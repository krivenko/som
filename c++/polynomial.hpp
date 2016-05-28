/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016 by I. Krivenko
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

#include <triqs/arrays/vector.hpp>

namespace som {

using triqs::arrays::vector;

class polynomial {

 // Coefficients of the polynomiall; the constant term goes first
 vector<double> coeffs;

public:

 polynomial(vector<double> const& coeffs = {}) : coeffs(coeffs) {}

 // Evaluation using Horner's rule
 double operator()(double x) const {
  double b = 0;
  for(int i = coeffs.size()-1; i>=0; --i) b = b*x + coeffs[i];
  return b;
 }

 double  operator[](int n) const { return coeffs[n]; }
 double& operator[](int n) { return coeffs[n]; }

};

}
