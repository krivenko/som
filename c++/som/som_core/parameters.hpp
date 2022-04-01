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

#include <cmath>
#include <functional>
#include <utility>
#include <variant>

#include <mpi/mpi.hpp>

#include <nda/nda.hpp>
#include <triqs/mesh/refreq.hpp>

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

  // Validate values of parameters
  void validate(observable_kind kind) const;
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

  /// Maximal ratio :math:`\chi/\chi_\mathrm{min}` for a particular solution to
  /// be considered good.
  double adjust_l_good_chi = 2.0;

  /// Maximal ratio :math:`\chi/\chi_\mathrm{min}` for a particular solution to
  /// be considered very good.
  double adjust_l_verygood_chi = 4.0 / 3.0;

  /// Critical ratio :math:`N_\mathrm{very good}/N_\mathrm{good}` to stop
  /// :math:`L`-adjustment procedure.
  double adjust_l_ratio = 0.95;

  /// Right boundary of the histograms, in units of :math:`\chi_\mathrm{min}`
  /// (left boundary is always set to :math:`\chi_\mathrm{min}`).
  double hist_max = 2.0;

  /// Number of bins for the histograms.
  int hist_n_bins = 100;

  accumulate_parameters_t() = default;
  explicit accumulate_parameters_t(std::pair<double, double> energy_window)
     : worker_parameters_t(energy_window) {}

  // Validate values of parameters
  void validate(observable_kind kind) const;
};

struct final_solution_cc_parameters_t {

  /// Grid of energy points used in derivative regularization procedure.
  std::variant<triqs::mesh::refreq, nda::array<double, 1>> refreq_mesh;

  /// Verbosity level (max level - 3).
  /// default: 2 on MPI rank 0, 0 otherwise.
  int verbosity =
      ((mpi::communicator().rank() == 0) ? 2 : 0); // silence the slave nodes

  /// Maximal ratio :math:`\chi/\chi_\mathrm{min}` for a particular solution to
  /// be selected. This criterion must be fulfilled together with the one set
  /// by `good_chi_abs`.
  double good_chi_rel = 2.0;

  /// Maximal value of :math:`\chi` for a particular solution to be selected.
  /// This criterion must be fulfilled together with the one set
  /// by `good_chi_rel`.
  double good_chi_abs = HUGE_VAL;

  /// Default model of the spectral function evaluated at
  /// energy points of `refreq_mesh`.
  nda::array<double, 1> default_model;

  /// Weights determining how much deviations from `default_model` are penalized
  /// at each energy point of `refreq_mesh`.
  nda::array<double, 1> default_model_weights;

  /// Maximum allowed number of parameter adjustment iterations.
  int max_iter = 20;

  /// Coefficient of the term that enforces the unity sum constraint.
  double unity_sum_coeff = 1e6;

  /// Maximum value of the regularization parameter that penalizes
  /// negative values of the spectral function.
  double amp_penalty_max = 1e3;

  /// Divisor used to reduce the regularization parameter that penalizes
  /// negative values of the spectral function.
  double amp_penalty_divisor = 10;

  /// Initial value of the regularization parameters that penalize large
  /// derivatives of the solution.
  double der_penalty_init = 1e-3;

  /// Coefficient used to increase the regularization parameters that penalize
  /// large derivatives of the solution.
  double der_penalty_coeff = 2.0;

  /// In the minimization of the CC quadratic form, a singular value
  /// :math:`\sigma_j` of its matrix is treated as zero if
  /// :math:`\sigma_j / \sigma_\mathrm{max}` is below this threshold.
  /// By default, the threshold is equal to the machine precision.
  double svd_rcond = -1;

  using monitor_t = std::function<bool(
      nda::vector<double>,
      std::pair<nda::array<double, 1>, nda::array<double, 1>>,
      std::pair<nda::array<double, 1>, nda::array<double, 1>>,
      std::pair<nda::array<double, 1>, nda::array<double, 1>>)>;

  /// Monitor function called at each parameter adjustment iteration.
  /// It takes the following arguments,
  /// - Current list of expansion coefficients :math:`c`;
  /// - Amplitudes of the spectral function and respective regularization
  ///   parameters as a :math:`(A_k, Q_k)` pair;
  /// - Derivatives of the spectral function and respective regularization
  ///   parameters as a :math:`(A'_k, D_k)` pair;
  /// - Second derivatives of the spectral function and respective
  ///   regularization parameters as a :math:`(A''_k, B_k)` pair.
  /// Returning `true` from the function stops the adjustment procedure.
  monitor_t monitor;

  final_solution_cc_parameters_t() = default;
  explicit final_solution_cc_parameters_t(triqs::mesh::refreq refreq_mesh)
     : refreq_mesh(std::move(refreq_mesh)) {}
  explicit final_solution_cc_parameters_t(nda::array<double, 1> refreq_mesh)
     : refreq_mesh(std::move(refreq_mesh)) {}

  // Validate values of parameters
  void validate() const;
};

} // namespace som
