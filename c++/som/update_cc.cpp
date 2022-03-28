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

#include <cmath>
#include <numeric>
#include <utility>

#ifdef EXT_DEBUG
#include <cassert>
#include <iostream>
#endif

#include <nda/linalg.hpp>

#include <triqs/utility/numeric_ops.hpp>

#include "numerics/ecqp_worker.hpp"
#include "numerics/finite_diff.hpp"

#include "configuration.hpp"
#include "kernels/all.hpp"
#include "solution_worker.hpp"
#include "update_cc.hpp"

namespace som {

using namespace triqs::arrays;

template <typename KernelType>
update_consistent_constraints<KernelType>::update_consistent_constraints(
    mc_data<KernelType>& data,
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
    double der_penalty_limiter)
   : data(data)
   , kern(data.objf.get_kernel())
   , rng(rng)
   , N(data.objf.get_rhs().size())
   , norm(norm)
   , norm_vector_view({1}, &(this->norm))
   , energy_window(std::move(energy_window))
   , window_width(energy_window.second - energy_window.first)
   , width_min(width_min)
   , weight_min(weight_min)
   , cycle_length(cycle_length)
   , max_rects(max_rects)
   , max_iter(max_iter)
   , rect_norm_variation_tol(rect_norm_variation_tol)
   , Q0_max(height_penalty_max * window_width / norm)
   , Q0_divisor(height_penalty_divisor)
   , Q1_init(der_penalty_init * std::pow(window_width, 2) / norm)
   , Q2_init(der_penalty_init * std::pow(window_width, 3) / norm)
   , Q12_threshold(der_penalty_threshold)
   , Q12_increase_coeff(der_penalty_increase_coeff)
   , Q12_limiter(der_penalty_limiter)
   , chi2_mat_prefactors(1.0 / (N * abs2(data.objf.get_error_bars())))
   , chi2_rhs_prefactors(chi2_mat_prefactors * conj(data.objf.get_rhs()))
#ifdef EXT_DEBUG
   , chi2_const(std::real(sum(chi2_rhs_prefactors * data.objf.get_rhs())))
#endif
   , proposed_conf(data.temp_conf.cache_ptr.get_ci()) {
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::operator()() {
  // CC-updates are proposed during the second stage of a global update
  // and only once in every cycle_length elementary updates.
  if(!(data.Z.in_stage_a() && (data.Z.current_step() % cycle_length == 0)))
    return;

  // Reject if a proposed non-overlapping configuration would be too big.
  // This check is conservative and can overestimate proposed_conf.size().
  if(2 * data.temp_conf.size() + 1 > max_rects) {
#ifdef EXT_DEBUG
    std::cerr
        << "Maximum number of rectangles in non-overlapping configuration "
        << (2 * data.temp_conf.size() + 1)
        << " would exceed max_rects = " << max_rects << ", rejecting"
        << std::endl;
#endif
    return;
  }

  ++n_proposed;
  double r = attempt();
  if(rng() < std::min(1.0, r)) {
    ++n_accepted;
    accept();
  } else
    reject();
}

template <typename KernelType>
double update_consistent_constraints<KernelType>::attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Proposing CC update" << std::endl;
#endif

  proposed_conf = make_nonoverlapping(data.temp_conf, energy_window, width_min);

  if(proposed_conf.size() < 3) {
#ifdef EXT_DEBUG
    std::cerr << proposed_conf.size()
              << " rectangles in the proposed configuration are too few "
                 "to apply the regularization functional, rejecting"
              << std::endl;
#endif
    return 0;
  }

#ifdef EXT_DEBUG
  std::cerr << "Optimizing heights of rectangles in non-overlapping "
               "configuration "
            << std::endl
            << proposed_conf << std::endl;
#endif

  auto K = long(proposed_conf.size());
  heights_prev.resize(K);
  heights.resize(K);
  proposed_conf.strip_rect_heights(heights);

  // Extract widths of the rectangles
  widths.resize(K);
  for(auto k : range(K)) { widths(k) = proposed_conf[k].width; }

  prepare_chi2_data(proposed_conf);

  optimize_heights();

  for(auto k : range(K)) {
    double rect_norm = heights(k) * widths(k);
    if(rect_norm < -weight_min) {
#ifdef EXT_DEBUG
      std::cerr << "CC optimization procedure produced a rectangle with a "
                   "significantly negative weight, rejecting (rectangle norm = "
                << rect_norm << ")" << std::endl;
#endif
      return 0;
    }
  }

  proposed_conf.update_rect_heights(heights);
  proposed_conf.redistribute_small_rects_weight(weight_min);

  new_objf_value = data.objf(proposed_conf);

#ifdef EXT_DEBUG
  std::cerr << "Final non-overlapping configuration: size = "
            << proposed_conf.size() << ", norm = " << norm
            << ", χ = " << std::sqrt(new_objf_value) << std::endl;
#endif

  return new_objf_value < data.temp_objf_value
             ? 1.0
             : data.Z(data.temp_objf_value / new_objf_value);
}

template <typename KernelType>
double update_consistent_constraints<KernelType>::accept() {
  // Update temporary configuration
  using std::swap;
  swap(proposed_conf, data.temp_conf);
  kern.cache_copy(proposed_conf, data.temp_conf);
  // and its objective function value
  data.temp_objf_value = new_objf_value;

#ifdef EXT_DEBUG
  std::cerr << "* CC update accepted" << std::endl;
  std::cerr << "Temporary configuration: size = " << data.temp_conf.size()
            << ", norm = " << data.temp_conf.norm()
            << ", χ = " << std::sqrt(data.temp_objf_value) << std::endl;

  for(auto& r : data.temp_conf) {
    if(r.norm() < weight_min)
      TRIQS_RUNTIME_ERROR << "Rectangle is too small: " << r;
    if(r.width < width_min)
      TRIQS_RUNTIME_ERROR << "Rectangle is too narrow: " << r;
    if(r.center - r.width / 2 < energy_window.first ||
       r.center + r.width / 2 > energy_window.second)
      TRIQS_RUNTIME_ERROR << "Rectangle is not within energy window: " << r;
  }
#endif

  // Copy this temporary configuration to the globally selected configuration
  // if its objective function value is smaller
  if(data.temp_objf_value < data.global_objf_value) {
#ifdef EXT_DEBUG
    std::cerr << "Copying temporary configuration to global configuration "
              << "(χ(temp) = " << std::sqrt(data.temp_objf_value)
              << ", χ(global) = " << std::sqrt(data.global_objf_value) << ")"
              << std::endl;
#endif
    data.global_conf = data.temp_conf;
    kern.cache_copy(data.temp_conf, data.global_conf);
    data.global_objf_value = data.temp_objf_value;
  }

#ifdef EXT_DEBUG
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

  // Update internal state of the distribution function
  ++data.Z;

  return 1;
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::reject() {

#ifdef EXT_DEBUG
  std::cerr << "* CC update rejected" << std::endl;
  std::cerr << "Temporary configuration: size = " << data.temp_conf.size()
            << ", norm = " << data.temp_conf.norm()
            << ", χ = " << std::sqrt(data.temp_objf_value) << std::endl;
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

  // Update internal state of the distribution function
  ++data.Z;
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::prepare_chi2_data(
    configuration const& conf) {
  auto K = long(conf.size());

  // Compute integrated kernels
  int_kernel.resize(K, N);
  for(auto k : range(K)) kern.apply(proposed_conf[k], int_kernel(k, range()));

  // Fill chi2_rhs
  chi2_rhs.resize(K);
  for(auto k : range(K)) {
    chi2_rhs(k) = std::real(sum(chi2_rhs_prefactors * int_kernel(k, range())));
  }

  // Fill chi2_mat
  chi2_mat.resize(K, K);
  for(auto k1 : range(K)) {
    for(auto k2 : range(K)) {
      chi2_mat(k1, k2) =
          std::real(sum(chi2_mat_prefactors * conj(int_kernel(k1, range())) *
                        int_kernel(k2, range())));
    }
  }

#ifdef EXT_DEBUG
  using triqs::utility::is_zero;
  if(!is_zero(max_element(abs(chi2_mat - transpose(chi2_mat))), 1e-10))
    TRIQS_RUNTIME_ERROR << "Matrix χ is not symmetric, χ = " << chi2_mat;
#endif
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::compute_O_mat() {
  auto K = long(proposed_conf.size());

  O_mat() = 0;

  // Contribution from the amplitude regularization functional
  // O_0 = sum_k Q_0(\epsilon_k)^2 |A(\epsilon_k)|^2
  for(auto k : range(K)) {
    double const Q0k = Q0(k);
    O_mat(k, k) += Q0k * Q0k;
  }

  // Contribution from the derivative regularization functional
  // O_1 = sum_k Q_1(\epsilon_k)^2 |A'(\epsilon_k)|^2
  for(auto k : range(K - 1)) {
    double const frac = Q1(k) / delta_c(k);
    double const val = frac * frac;
    O_mat(k, k) += val;
    O_mat(k + 1, k + 1) += val;
    O_mat(k, k + 1) += -val;
    O_mat(k + 1, k) += -val;
  }

  // Contribution from the second derivative regularization functional
  // O_2 = sum_k Q_2(\epsilon_k)^2 |A''(\epsilon_k)|^2
  for(auto k : range(1, K - 1)) {
    double const Q2k = Q2(k - 1);
    double const coeff = 4 * Q2k * Q2k;
    double const delta_km1 = delta_c(k - 1);
    double const delta_kp1 = delta_c(k);
    double const delta_k = delta_km1 + delta_kp1;
    O_mat(k - 1, k - 1) += coeff / (delta_km1 * delta_km1 * delta_k * delta_k);
    O_mat(k + 1, k + 1) += coeff / (delta_kp1 * delta_kp1 * delta_k * delta_k);
    O_mat(k, k) += coeff / (delta_km1 * delta_km1 * delta_kp1 * delta_kp1);
    O_mat(k - 1, k + 1) += coeff / (delta_km1 * delta_kp1 * delta_k * delta_k);
    O_mat(k + 1, k - 1) += coeff / (delta_km1 * delta_kp1 * delta_k * delta_k);
    O_mat(k - 1, k) += -coeff / (delta_km1 * delta_km1 * delta_kp1 * delta_k);
    O_mat(k, k - 1) += -coeff / (delta_km1 * delta_km1 * delta_kp1 * delta_k);
    O_mat(k, k + 1) += -coeff / (delta_km1 * delta_kp1 * delta_kp1 * delta_k);
    O_mat(k + 1, k) += -coeff / (delta_km1 * delta_kp1 * delta_kp1 * delta_k);
  }

#ifdef EXT_DEBUG
  using triqs::utility::is_zero;
  if(!is_zero(max_element(abs(O_mat - transpose(O_mat))), 1e-10))
    TRIQS_RUNTIME_ERROR << "Matrix O is not symmetric, O = " << O_mat;
#endif
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::optimize_heights() {
  auto K = long(proposed_conf.size());

  // Compute distances between centers of adjacent rectangles
  delta_c.resize(K - 1);
  for(auto k : range(K - 1)) {
    delta_c(k) = proposed_conf[k + 1].center - proposed_conf[k].center;
  }

  Q0.resize(K);
  Q0() = 0;
  Q1.resize(K - 1);
  Q1() = Q1_init;
  Q2.resize(K - 2);
  Q2() = Q2_init;

  O_mat.resize(K, K);
  QF_mat.resize(K, K);

  conf_1st_der.resize(K - 1);
  conf_2nd_der.resize(K - 2);

  double const C = Q12_threshold / std::sqrt(K - 1);

  auto worker = ecqp_worker(int(K), 1);
  // View of the 1D array 'widths' as a 1xK matrix.
  auto widths_matrix_view =
      nda::matrix_const_view<double, nda::F_layout>({1, K}, widths.data());

  int iter = 0;
  for(; iter < max_iter; ++iter) {

    std::swap(heights_prev, heights);

    // Update the heights
    compute_O_mat();
    QF_mat() = chi2_mat + O_mat;
#ifdef EXT_DEBUG
    double QF =
#endif
        worker(QF_mat, chi2_rhs, widths_matrix_view, norm_vector_view, heights);

    double norm_diff_max =
        sum(abs((make_array_view(heights) - make_array_view(heights_prev)) *
                widths)) /
        norm;

#ifdef EXT_DEBUG
    double O = nda::blas::dot(heights, O_mat * heights);
    double chi2 = QF + chi2_const - O;

    std::cerr << "  Iteration " << (iter + 1) << " out of " << max_iter << ": "
              << "sum(|Δnorm_k|) / norm = " << norm_diff_max
              << ", χ^2 = " << chi2 << ", O/χ^2 = " << (O / chi2) << std::endl;
#endif

    if(norm_diff_max < rect_norm_variation_tol) {
#ifdef EXT_DEBUG
      nda::array<double, 1> rectangle_norms = make_array_view(heights) * widths;
      std::cerr << "Convergence reached, rectangle norms = " << rectangle_norms
                << std::endl;
#endif
      break;
    }

    // Update amplitude regularization parameters
    for(auto k : range(K))
      Q0(k) = (heights(k) < 0) ? Q0_max : (Q0(k) / Q0_divisor);

    // Update 1st derivative regularization parameters
    finite_diff_forward_dx(heights(), delta_c(), conf_1st_der());

    double Q1_limit = Q12_limiter * min_element(Q1);
    for(auto k : range(K - 1)) {
      double d = std::abs(conf_1st_der(k));
      if(d > C / Q1(k))
        Q1(k) = C / d;
      else
        Q1(k) *= Q12_increase_coeff;
      if(Q1(k) > Q1_limit) Q1(k) = Q1_limit;
    }

    // Update 2nd derivative regularization parameters
    finite_diff_2_symm_dx(heights(), delta_c(), conf_2nd_der());

    double Q2_limit = Q12_limiter * min_element(Q2);
    for(auto k : range(K - 2)) {
      double d = std::abs(conf_2nd_der(k));
      if(d > C / Q2(k))
        Q2(k) = C / d;
      else
        Q2(k) *= Q12_increase_coeff;
      if(Q2(k) > Q2_limit) Q2(k) = Q2_limit;
    }
  }

#ifdef EXT_DEBUG
  if(iter == max_iter) {
    std::cerr << "Maximum number of iterations reached, heights = " << heights
              << std::endl;
  }
#endif
}

template <typename KernelType>
double update_consistent_constraints<KernelType>::get_acceptance_rate() const {
  return double(n_accepted) / double(n_proposed);
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_consistent_constraints);

} // namespace som
