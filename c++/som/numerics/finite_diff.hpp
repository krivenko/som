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
#pragma once

#include <cassert>

#include <triqs/arrays.hpp>

// Computation of finite differences

namespace som {

// Forward difference
template <typename T1, typename T2>
inline void finite_diff_forward(triqs::arrays::array_const_view<T1, 1> f,
                                triqs::arrays::array<double, 1> const& x,
                                triqs::arrays::array_view<T2, 1> output) {
  assert(f.size() == x.size());
  assert(output.size() == x.size() - 1);
  for(int k : triqs::arrays::range(output.size())) {
    output(k) = (f(k + 1) - f(k)) / (x(k + 1) - x(k));
  }
}

// Second order symmetric difference
template <typename T1, typename T2>
inline void finite_diff_2_symm(triqs::arrays::array_const_view<T1, 1> f,
                               triqs::arrays::array<double, 1> const& x,
                               triqs::arrays::array_view<T2, 1> output) {
  assert(f.size() == x.size());
  assert(output.size() == x.size() - 2);
  for(int k : triqs::arrays::range(output.size())) {
    output(k) =
        (f(k + 2) * (x(k + 1) - x(k)) + f(k) * (x(k + 2) - x(k + 1)) -
         f(k + 1) * (x(k + 2) - x(k))) /
        (0.5 * (x(k + 2) - x(k)) * (x(k + 2) - x(k + 1)) * (x(k + 1) - x(k)));
  }
}

} // namespace som
