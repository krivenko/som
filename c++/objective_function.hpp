/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2017 by I. Krivenko
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
#include "config_update.hpp"

namespace som {

template<typename KernelType>
class objective_function {

 using rhs_type = typename KernelType::result_type;
 using mesh_type = typename KernelType::mesh_type;

 // Integral kernel
 KernelType const& kern;
 // The right-hand side of the Fredholm integral equation
 rhs_type const& rhs;
 // Error bars of the RHS
 rhs_type const& error_bars;

public:

 objective_function(KernelType const& kern,
                    rhs_type const& rhs,
                    rhs_type const& error_bars) :
  kern(kern), rhs(rhs), error_bars(error_bars) {}

 double operator()(configuration const& c) const {
  return sum(abs((rhs - kern(c)) / error_bars));
 }

 double operator()(config_update const& cu) const {
  return sum(abs((rhs - kern(cu)) / error_bars));
 }

 KernelType const& get_kernel() const { return kern; }
};

}
