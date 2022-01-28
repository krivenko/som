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

#include <utility>

#include <som/worker_parameters.hpp>

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

// Arguments of the som_core::accumulate() function
struct accumulate_parameters_t : public worker_parameters_t {

  /////////////////////
  // Main parameters //
  /////////////////////

  /// Number of global updates (:math:`F`);
  int f = 100;

  /// Number of particular solutions used in accumulation (:math:`L`);
  /// Ignored if `adjust_l = True`.
  int l = 2000;

  /// Adjust the number of solutions used in accumulation.
  bool adjust_l = false;

  /// Accumulate histograms of objective function values.
  bool make_histograms = false;

  /// Search range for the number of solutions used in the final accumulation.
  /// type: (int,int)
  std::pair<int, int> adjust_l_range = std::pair<int, int>{100, 2000};

  /// Maximal ratio :math:`D/D_\mathrm{min}` for a particular solution to be
  /// considered good.
  double adjust_l_good_d = 2.0;

  /// Maximal ratio :math:`D/D_\mathrm{min}` for a particular solution to be
  /// considered very good.
  double adjust_l_verygood_d = 4.0 / 3.0;

  /// Critical ratio :math:`N_\mathrm{very good}/N_\mathrm{good}` to stop
  /// :math:`L`-adjustment procedure.
  double adjust_l_ratio = 0.95;

  /// Right boundary of the histograms, in units of :math:`D_\mathrm{min}`
  /// (left boundary is always set to :math:`D_\mathrm{min}`).
  double hist_max = 2.0;

  /// Number of bins for the histograms.
  int hist_n_bins = 100;

  accumulate_parameters_t() = default;
  explicit accumulate_parameters_t(std::pair<double, double> energy_window)
     : worker_parameters_t(energy_window) {}
};

} // namespace som