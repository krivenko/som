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
#pragma once

#include <cmath>

namespace som {

// This function object returns a real random number x from [x_min; x_max[
// distributed as P(x) = N exp(-gamma*|x|/max(|x_min|,|x_max|))
template <typename RNG> class expabs_distribution {

  RNG& rng;
  const double gamma;

public:
  expabs_distribution(RNG& rng, double gamma) : rng(rng), gamma(gamma) {}

  double operator()(double x_min, double x_max) const {
    double x_min_abs = std::abs(x_min);
    double x_max_abs = std::abs(x_max);
    double L = std::max(x_min_abs, x_max_abs);
    double gL = gamma / L;

    double x = rng();

    double f = (1 - x) * std::copysign(std::expm1(-gL * x_min_abs), x_min) +
               x * std::copysign(std::expm1(-gL * x_max_abs), x_max);

    return std::copysign(std::log1p(-std::abs(f)) / gL, f);
  }
};

} // namespace som
