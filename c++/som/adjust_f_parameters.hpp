/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <utility>

#include "worker_parameters.hpp"

namespace som {

// Arguments of the som_core::adjust_f() function
struct adjust_f_parameters_t : public worker_parameters_t {

  /// Search range for the number of global updates.
  /// type: (int,int)
  std::pair<int, int> f_range = std::pair<int, int>{100, 5000};

  /// Number of particular solutions used to adjust :math:`F`.
  int l = 20;

  /// Limiting value of :math:`\kappa` used to adjust :math:`F`.
  double kappa = 0.25;

  adjust_f_parameters_t() = default;
  explicit adjust_f_parameters_t(std::pair<double, double> energy_window)
     : worker_parameters_t(energy_window) {}
};

} // namespace som
