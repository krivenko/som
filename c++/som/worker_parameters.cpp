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

#include <string>
#include <utility>

#include <triqs/utility/exceptions.hpp>

#include <som/kernels/observables.hpp>

#include "worker_parameters.hpp"

namespace som {

[[noreturn]] void fatal_error(std::string const& message) {
  TRIQS_RUNTIME_ERROR << "solution_worker: " << message;
}

void worker_parameters_t::validate(observable_kind kind) const {
  using std::to_string;

  if(energy_window.first >= energy_window.second)
    fatal_error("Wrong energy window [" + to_string(energy_window.first) + ";" +
                to_string(energy_window.second) + "]");

  auto [e_min, e_max] = max_energy_window(kind);
  if(energy_window.first < e_min) {
    fatal_error("Left boundary of the energy window " +
                to_string(energy_window.first) +
                " is invalid for this observable");
  }
  if(energy_window.second > e_max) {
    fatal_error("Right boundary of the energy window " +
                to_string(energy_window.second) +
                " is invalid for this observable");
  }

  if(t <= 0)
    fatal_error("Number of elementary updates t must be positive (got t = " +
                to_string(t) + ")");

  if(max_rects <= 0)
    fatal_error(
        "Maximum number of rectangles must be positive (got max_rects = " +
        to_string(max_rects) + ")");

  if(min_rect_width <= 0 || min_rect_width >= 1)
    fatal_error("min_rect_width = " + to_string(min_rect_width) +
                " must be in (0;1)");

  if(min_rect_weight <= 0 || min_rect_weight > 1)
    fatal_error("min_rect_weight = " + to_string(min_rect_weight) +
                " must be in (0;1]");

  if(distrib_d_max < 1)
    fatal_error("Parameter distrib_d_max = " + to_string(distrib_d_max) +
                " must be at least 1");

  if(gamma <= 0)
    fatal_error(
        "Proposal probability parameter gamma must be positive"
        " (got gamma = " +
        to_string(gamma) + ")");
}

} // namespace som
