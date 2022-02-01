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

#include <som/config_update.hpp>
#include <som/configuration.hpp>
#include <som/kernels/all.hpp>

namespace som {

template <typename KernelType> class objective_function {

  using rhs_type = typename KernelType::result_type;
  using mesh_type = typename KernelType::mesh_type;

  // Integral kernel
  KernelType const& kern;
  // The right-hand side of the Fredholm integral equation
  rhs_type const& rhs;
  // Error bars of the RHS
  rhs_type const& error_bars;

public:
  objective_function(KernelType const& kern, rhs_type const& rhs,
                     rhs_type const& error_bars);

  double operator()(configuration const& c) const;
  double operator()(config_update const& cu) const;

  [[nodiscard]] KernelType const& get_kernel() const { return kern; }
};

EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(objective_function)

} // namespace som
