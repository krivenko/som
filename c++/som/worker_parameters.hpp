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

  /// Maximal parameter of the power-law distribution function for the
  /// Metropolis algorithm.
  double distrib_d_max = 2;

  /// Proposal probability parameter :math:`\gamma`.
  double gamma = 2;

  worker_parameters_t() = default;
  explicit worker_parameters_t(std::pair<double, double> energy_window)
     : energy_window(std::move(energy_window)) {}
};

} // namespace som
