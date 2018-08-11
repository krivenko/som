/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2018 by I. Krivenko
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

#include <vector>
#include <cmath>
#include <utility>

#include <triqs/arrays/vector.hpp>

#include "configuration.hpp"

namespace som {

// Fit quality function (kappa)
template<typename KernelType>
class fit_quality {

 using rhs_type = typename KernelType::result_type;
 using mesh_type = typename KernelType::mesh_type;

 // Integral kernel
 KernelType const& kern;
 // The right-hand side of the Fredholm integral equation
 rhs_type const& rhs;
 // Error bars of the RHS
 rhs_type const& error_bars;

 // Correlation measure of two real numbers
 static double corr(double x1, double x2) {
  return (x1*x2 < 0) ? 1 : 0;
 }

 // Correlation measure of two complex numbers
 static double corr(dcomplex z1, dcomplex z2) {
  dcomplex w = z1 * std::conj(z2);
  return 0.5*(1.0 - (w == .0 ? 1.0 : (w.real() / std::abs(w))));
 }

public:

 fit_quality(KernelType const& kern,
             rhs_type const& rhs,
             rhs_type const& error_bars) :
  kern(kern), rhs(rhs), error_bars(error_bars) {}

 double operator()(configuration const& c) const {
  auto delta = (rhs - kern(c)) / error_bars;
  int M = first_dim(delta);
  double kappa = 0;
  for(int i = 1; i < M; ++i) kappa += corr(delta(i), delta(i-1));
  kappa /= M - 1;
  return kappa;
 }

 KernelType const& get_kernel() const { return kern; }
};

}
