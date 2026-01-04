/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2026 Igor Krivenko
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

[[noreturn]] inline static void fatal_error(std::string const& message) {
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

  if(use_symmetrized_spectrum(kind) &&
     (energy_window.first != -energy_window.second))
    fatal_error(
        "The energy window must be symmetric w.r.t. zero "
        "for this observable");

  if(t <= 0)
    fatal_error("Number of elementary updates t must be positive (got t = " +
                to_string(t) + ")");

  if(max_rects <= 0)
    fatal_error(
        "Maximum number of rectangles must be positive (got max_rects = " +
        to_string(max_rects) + ")");

  if(min_rect_width < 0 || min_rect_width >= 1)
    fatal_error("min_rect_width = " + to_string(min_rect_width) +
                " must be in [0;1[");

  if(min_rect_weight < 0 || min_rect_weight > 1)
    fatal_error("min_rect_weight = " + to_string(min_rect_weight) +
                " must be in [0;1]");

  if(t1 < -1 || t1 > t)
    fatal_error("Number of updates in the first stage, t1 = " + to_string(t1) +
                " must be in [0;t] = [0;" + to_string(t) + "]");

  if(distrib_d_max < 1)
    fatal_error("Parameter distrib_d_max = " + to_string(distrib_d_max) +
                " must be at least 1");

  if(gamma <= 0)
    fatal_error(
        "Proposal probability parameter gamma must be positive"
        " (got gamma = " +
        to_string(gamma) + ")");

  if(cc_update) {
    if(cc_update_cycle_length <= 0)
      fatal_error("CC parameter cc_update_cycle_length must be positive, got " +
                  to_string(cc_update_cycle_length));

    if(cc_update_max_iter <= 0)
      fatal_error("CC parameter cc_update_max_iter must be positive, got " +
                  to_string(cc_update_max_iter));

    if(cc_update_rect_norm_variation_tol <= 0)
      fatal_error(
          "CC parameter cc_update_rect_norm_variation_tol must be positive, "
          "got " +
          to_string(cc_update_rect_norm_variation_tol));

    if(cc_update_height_penalty_max <= 0)
      fatal_error(
          "CC parameter cc_update_height_penalty_max must be positive, got " +
          to_string(cc_update_height_penalty_max));

    if(cc_update_height_penalty_divisor <= 1)
      fatal_error(
          "CC parameter cc_update_height_penalty_divisor must exceed 1.0, "
          "got " +
          to_string(cc_update_height_penalty_divisor));

    if(cc_update_der_penalty_init <= 0)
      fatal_error(
          "CC parameter cc_update_der_penalty_init must be positive, got " +
          to_string(cc_update_der_penalty_init));

    if(cc_update_der_penalty_threshold <= 0)
      fatal_error(
          "CC parameter cc_update_der_penalty_threshold must be positive, "
          "got " +
          to_string(cc_update_der_penalty_threshold));

    if(cc_update_der_penalty_increase_coeff <= 1)
      fatal_error(
          "CC parameter cc_update_der_penalty_increase_coeff must exceed 1.0, "
          "got " +
          to_string(cc_update_der_penalty_increase_coeff));

    if(cc_update_der_penalty_limiter <= 1)
      fatal_error(
          "CC parameter cc_update_der_penalty_limiter must exceed 1.0, got " +
          to_string(cc_update_der_penalty_limiter));

    if(cc_update_rect_norm_variation_tol > min_rect_weight)
      std::cout << "WARNING: cc_update_rect_norm_variation_tol should normally "
                   "be set below min_rect_weight\n";
  }
}

} // namespace som
