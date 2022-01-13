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

#include <triqs/arrays/vector.hpp>

#include "objective_function.hpp"

namespace som {

template <typename KernelType>
objective_function<KernelType>::objective_function(KernelType const& kern,
                                                   rhs_type const& rhs,
                                                   rhs_type const& error_bars)
   : kern(kern), rhs(rhs), error_bars(error_bars) {}

template <typename KernelType>
double
objective_function<KernelType>::operator()(configuration const& c) const {
  return sum(abs((rhs - kern(c)) / error_bars));
}

template <typename KernelType>
double
objective_function<KernelType>::operator()(config_update const& cu) const {
  return sum(abs((rhs - kern(cu)) / error_bars));
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(objective_function)

} // namespace som
