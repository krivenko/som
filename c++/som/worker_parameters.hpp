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

#include <string>
#include <utility>

#include <mpi/mpi.hpp>

namespace som {

enum observable_kind : unsigned int;

/// Parameters of solution_worker::operator()()
struct worker_parameters_t {

  /////////////////////
  // Main parameters //
  /////////////////////

  /// Estimated lower and upper bounds of the spectrum.
  /// Negative values of the lower bound will be reset to 0 for susceptibilities
  /// and conductivity.
  /// type: (float,float)
  std::pair<double, double> energy_window;

  /// Maximum runtime in seconds, use -1 to set infinite.
  /// default: -1 = infinite
  int max_time = -1;

  /// Verbosity level (max level - 3).
  /// default: 2 on MPI rank 0, 0 otherwise.
  int verbosity =
      ((mpi::communicator().rank() == 0) ? 2 : 0); // silence the slave nodes

  /// Number of elementary updates per global update (:math:`T`).
  int t = 50;

  /// Enable Consistent Constraints updates.
  bool cc_update = false;

  /////////////////////////
  // Fine tuning options //
  /////////////////////////

  /// Seed for random number generator.
  /// default: 34788 + 928374 * MPI.rank
  int random_seed = 34788 + 928374 * mpi::communicator().rank();

  /// Name of random number generator.
  /// type: str
  std::string random_name = "";

  /// Maximum number of rectangles to represent spectra (:math:`K_{max}`),
  /// should be below 70.
  int max_rects = 60;

  /// Minimal width of a rectangle, in units of the energy window width.
  double min_rect_width = 1e-3;

  /// Minimal weight of a rectangle, in units of the requested solution norm.
  double min_rect_weight = 1e-3;

  /// Number of elementary updates in the first stage of a global update.
  /// When set to -1 (default), the number of elementary updates will be chosen
  /// randomly for each global update from the :math:`[1; T[` range.
  int t1 = -1;

  /// Maximal parameter of the power-law distribution function for the
  /// Metropolis algorithm.
  double distrib_d_max = 2;

  /// Proposal probability parameter :math:`\gamma`.
  double gamma = 2;

  /// CC update: Number of proposed elementary updates between two successive
  /// CC updates (only during stage A of a global update).
  int cc_update_cycle_length = 10;

  /// CC update: Maximum allowed number of height adjustment iterations.
  int cc_update_max_iter = 30;

  /// CC update: The height adjustment procedure stops when variation of every
  /// rectangle norm between two consecutive iterations is below this value.
  /// This parameter is measured in units of the requested solution norm.
  double cc_update_rect_norm_variation_tol = 1e-3;

  /// CC update: Maximum value of the regularization parameters :math:`Q_0(k)`
  /// that penalize negative heights.
  /// Measured in units of (energy window width) / (solution norm).
  double cc_update_height_penalty_max = 1e3;

  /// CC update: Divisor used to reduce the regularization parameters
  /// :math:`Q_0(k)` that penalize negative heights.
  double cc_update_height_penalty_divisor = 10;

  /// CC update: Initial value of the regularization parameters :math:`Q_1(k)`
  /// and :math:`Q_2(k)` that penalize large derivatives of a solution.
  /// Measured in units of (energy window width)^2 / (solution norm) for
  /// :math:`Q_1(k)` and in units of (energy window width)^3 / (solution norm)
  /// for :math:`Q_2(k)`.
  double cc_update_der_penalty_init = 1;

  /// CC update: Sets the threshold value of the products :math:`|Q_1(k) A'(k)|`
  /// and :math:`Q_2(k) A''(k)`, above which derivative regularization
  /// parameters :math:`Q_1(k)` and :math:`Q_2(k)` need to be reduced.
  double cc_update_der_penalty_threshold = 0.1;

  /// CC update: Coefficient used to increase the regularization parameters
  /// :math:`Q_1(k)` and :math:`Q_2(k)` that penalize large derivatives of
  /// a solution.
  double cc_update_der_penalty_increase_coeff = 2;

  /// CC update: Coefficient that limits growth of the regularization parameters
  /// :math:`Q_1(k)` and :math:`Q_2(k)` that penalize large derivatives of
  /// a solution.
  double cc_update_der_penalty_limiter = 1e3;

  worker_parameters_t() = default;
  explicit worker_parameters_t(std::pair<double, double> energy_window)
     : energy_window(std::move(energy_window)) {}

  // Validate values of parameters
  void validate(observable_kind kind) const;
};

} // namespace som
