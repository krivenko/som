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

#include <nda/nda.hpp>

#include <triqs/mc_tools/random_generator.hpp>

#include "solution_functionals/objective_function.hpp"

namespace som {

template <typename> struct mc_data;

// Consistent Constraints update
template <typename KernelType> class update_consistent_constraints {

  using rhs_type = typename KernelType::result_type;

  using rhs_scalar_type = typename rhs_type::value_type;

  mc_data<KernelType>& data;
  KernelType const& kern;
  triqs::mc_tools::random_generator& rng;

  int const N_sigma2;
  double const norm;
  nda::vector_const_view<double> const norm_vector_view;

  std::pair<double, double> const energy_window;
  double const window_width;
  double const width_min;
  double const weight_min;

  int const cycle_length;
  std::uint64_t n_proposed = 0;
  std::uint64_t n_accepted = 0;

  int const max_rects;

  int const max_iter;
  double const rect_norm_variation_tol;
  double const Q0_max;
  double const Q0_divisor;
  double const Q1_init;
  double const Q2_init;
  double const Q12_threshold;
  double const Q12_increase_coeff;
  double const Q12_limiter;

  // Coefficients (1 / \tilde N) (\sigma_n^2 + l^2)^{-1}
  nda::array<double, 1> const chi2_eigenvalue_prefactors;
  // conj(\hat U^\dagger g_n)
  nda::array<rhs_scalar_type, 1> const chi2_conj_rhs;
#ifdef EXT_DEBUG
  // Constant (1/\tilde N) \sum_n |\hat U^\dagger g_n|^2 / (\sigma_n^2 + l^2)
  double const chi2_const;
#endif

  // Proposed non-overlapping configuration.
  configuration proposed_conf;

  // range(K), range(K - 1), range(K - 2)
  nda::range r_K, r_K1, r_K2;

  // Objective function value for the proposed configuration.
  double new_objf_value = NAN;

  // Kernel, integrated over one unity-rectangle.
  nda::array<rhs_scalar_type, 1> int_kernel_one_rect;

  // U^\dagger-transformed kernel, integrated over all unity-rectangles in
  // the proposed configuration.#include <iomanip> //FIXME
  nda::array<rhs_scalar_type, 2> int_kernel;

  // Matrix of quadratic form \chi^2.
  nda::matrix<double, nda::F_layout> chi2_mat;
  // Linear shifts in \chi^2.
  nda::vector<double> chi2_rhs;

  // Widths of rectangles
  nda::array<double, 1> widths;

  // Distances between centers of adjacent rectangles.
  nda::array<double, 1> delta_c;

  // Heights of rectangles on the previous and current CC optimization
  // iteration.
  nda::vector<double> heights_prev, heights;

  // Parameters of amplitude regularization functional O_0.
  nda::array<double, 1> Q0;
  // Parameters of derivative regularization functional O_1.
  nda::array<double, 1> Q1;
  // Parameters of second derivative regularization functional O_2.
  nda::array<double, 1> Q2;

  // Matrix of regularizing quadratic form O.
  nda::matrix<double, nda::F_layout> O_mat;

  // Matrix of quadratic form \chi^2 + O.
  // This array serves as a preallocated buffer of size max_rects * max_rects,
  // matrix view of it with shape (K, K) is constructed in optimize_heights().
  nda::array<double, 1> QF_mat_data;

  // First derivative of the proposed configuration.
  nda::array<double, 1> conf_1st_der;
  // Second derivative of the proposed configuration.
  nda::array<double, 1> conf_2nd_der;

  // Initialize chi2_conj_rhs
  static nda::array<rhs_scalar_type, 1>
  init_chi2_conj_rhs(objective_function<KernelType> const& objf);

  // Prepare matrix of quadratic form \chi^2 and coefficients of
  // linear shifts in \chi^2.
  void prepare_chi2_data();

  // Compute matrix of regularizing quadratic form O.
  void compute_O_mat();

  // Self-consistent optimization of rectangle heights.
  void optimize_heights();

public:
  update_consistent_constraints(mc_data<KernelType>& data,
                                triqs::mc_tools::random_generator& rng,
                                double norm,
                                std::pair<double, double> energy_window,
                                double width_min,
                                double weight_min,
                                int cycle_length,
                                int max_rects,
                                int max_iter,
                                double rect_norm_variation_tol,
                                double height_penalty_max,
                                double height_penalty_divisor,
                                double der_penalty_init,
                                double der_penalty_threshold,
                                double der_penalty_increase_coeff,
                                double der_penalty_limiter);

  // Call attempt() and accept()/reject() according to the Metropolis criterion
  void operator()();

  double attempt();
  double accept();
  void reject();

  std::uint64_t get_n_proposed() const { return n_proposed; }
  double get_acceptance_rate() const;
};

EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(update_consistent_constraints);

} // namespace som
