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
   , N_sigma2(data.objf.get_sigma2().size())
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
   , chi2_eigenvalue_prefactors(1.0 / (N_sigma2 * data.objf.get_sigma2()))
   , chi2_conj_rhs(init_chi2_conj_rhs(data.objf))
#ifdef EXT_DEBUG
   , chi2_const(
         std::real(sum(chi2_eigenvalue_prefactors * abs2(chi2_conj_rhs))))
#endif
   , proposed_conf(data.temp_conf.cache_ptr.get_ci())
   , int_kernel_one_rect(data.objf.get_rhs().size())
   , int_kernel(max_rects, N_sigma2)
   , chi2_mat(max_rects, max_rects)
   , chi2_rhs(max_rects)
   , widths(max_rects)
   , delta_c(max_rects)
   , heights_prev(max_rects)
   , heights(max_rects)
   , Q0(max_rects)
   , Q1(max_rects - 1)
   , Q2(max_rects - 2)
   , O_mat(max_rects, max_rects)
   , QF_mat_data(max_rects * max_rects)
   , conf_1st_der(max_rects - 1)
   , conf_2nd_der(max_rects - 2) {
}

template <typename KernelType>
auto update_consistent_constraints<KernelType>::init_chi2_conj_rhs(
    objective_function<KernelType> const& objf)
    -> nda::array<rhs_scalar_type, 1> {
  auto const& U_dagger = objf.get_U_dagger();
  if(U_dagger)
    return conj((*U_dagger) *
                vector_const_view<rhs_scalar_type>(objf.get_rhs()));
  else
    return conj(objf.get_rhs());
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
        << " would exceed max_rects = " << max_rects << ", rejecting\n";
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
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  std::cerr << "* Proposing CC update\n";
#endif

  proposed_conf = make_nonoverlapping(data.temp_conf, energy_window, width_min);

  if(proposed_conf.size() < 3) {
#ifdef EXT_DEBUG
    std::cerr << proposed_conf.size()
              << " rectangles in the proposed configuration are too few "
                 "to apply the regularization functional, rejecting\n";
#endif
    return 0;
  }

#ifdef EXT_DEBUG
  std::cerr << "Optimizing heights of rectangles in non-overlapping "
               "configuration \n"
            << proposed_conf << '\n';
#endif

  auto K = long(proposed_conf.size());
  r_K = range(K);
  r_K1 = range(K - 1);
  r_K2 = range(K - 2);

  proposed_conf.strip_rect_heights(heights(r_K));

  // Extract widths of the rectangles
  for(auto k : r_K) { widths(k) = proposed_conf[k].width; }

  prepare_chi2_data();

  optimize_heights();

  for(auto k : r_K) {
    double rect_norm = heights(k) * widths(k);
    if(rect_norm < -weight_min) {
#ifdef EXT_DEBUG
      std::cerr << "CC optimization procedure produced a rectangle with a "
                   "significantly negative weight, rejecting (rectangle norm = "
                << rect_norm << ")\n";
#endif
      return 0;
    }
  }

  proposed_conf.update_rect_heights(heights(r_K));
  proposed_conf.redistribute_small_rects_weight(weight_min);

  new_objf_value = data.objf(proposed_conf);

#ifdef EXT_DEBUG
  std::cerr << "Final non-overlapping configuration: size = "
            << proposed_conf.size() << ", norm = " << norm
            << ", χ = " << std::sqrt(new_objf_value) << '\n';
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
  std::cerr << "* CC update accepted\n";
  std::cerr << "Temporary configuration: size = " << data.temp_conf.size()
            << ", norm = " << data.temp_conf.norm()
            << ", χ = " << std::sqrt(data.temp_objf_value) << '\n';

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
              << ", χ(global) = " << std::sqrt(data.global_objf_value) << ")\n";
#endif
    data.global_conf = data.temp_conf;
    kern.cache_copy(data.temp_conf, data.global_conf);
    data.global_objf_value = data.temp_objf_value;
  }

#ifdef EXT_DEBUG
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
#endif

  // Update internal state of the distribution function
  ++data.Z;

  return 1;
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::reject() {

#ifdef EXT_DEBUG
  std::cerr << "* CC update rejected\n";
  std::cerr << "Temporary configuration: size = " << data.temp_conf.size()
            << ", norm = " << data.temp_conf.norm()
            << ", χ = " << std::sqrt(data.temp_objf_value) << '\n';
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
#endif

  // Update internal state of the distribution function
  ++data.Z;
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::prepare_chi2_data() {
  // Compute integrated kernels
  auto const& U_dagger = data.objf.get_U_dagger();
  if(U_dagger) { // Covariance matrix
    for(auto k : r_K) {
      kern.apply(proposed_conf[k], int_kernel_one_rect());
      int_kernel(k, range::all) =
          (*U_dagger) * vector_const_view<rhs_scalar_type>(int_kernel_one_rect);
    }
  } else { // Estimated error bars
    for(auto k : r_K) kern.apply(proposed_conf[k], int_kernel(k, range::all));
  }

  // Fill chi2_rhs
  for(auto k : r_K) {
    chi2_rhs(k) = std::real(sum(chi2_eigenvalue_prefactors *
                                int_kernel(k, range::all) * chi2_conj_rhs));
  }

  // Fill chi2_mat
  for(auto k1 : r_K) {
    for(auto k2 : r_K) {
      chi2_mat(k1, k2) = std::real(
          sum(chi2_eigenvalue_prefactors * conj(int_kernel(k1, range::all)) *
              int_kernel(k2, range::all)));
    }
  }

#ifdef EXT_DEBUG
  using triqs::utility::is_zero;
  auto chi2_mat_view = chi2_mat(r_K, r_K);
  if(!is_zero(max_element(abs(chi2_mat_view - transpose(chi2_mat_view))),
              1e-10))
    TRIQS_RUNTIME_ERROR << "Matrix χ is not symmetric, χ = " << chi2_mat_view;
#endif
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::compute_O_mat() {
  O_mat(r_K, r_K) = 0;

  // Contribution from the amplitude regularization functional
  // O_0 = sum_k Q_0(\epsilon_k)^2 |A(\epsilon_k)|^2
  for(auto k : r_K) {
    double const Q0k = Q0(k);
    O_mat(k, k) += Q0k * Q0k;
  }

  // Contribution from the derivative regularization functional
  // O_1 = sum_k Q_1(\epsilon_k)^2 |A'(\epsilon_k)|^2
  for(auto k : r_K1) {
    double const frac = Q1(k) / delta_c(k);
    double const val = frac * frac;
    O_mat(k, k) += val;
    O_mat(k + 1, k + 1) += val;
    O_mat(k, k + 1) += -val;
    O_mat(k + 1, k) += -val;
  }

  // Contribution from the second derivative regularization functional
  // O_2 = sum_k Q_2(\epsilon_k)^2 |A''(\epsilon_k)|^2
  for(auto k : r_K2) {
    double const Q2k = Q2(k);
    double const coeff = 4 * Q2k * Q2k;
    double const delta_km1 = delta_c(k);
    double const delta_kp1 = delta_c(k + 1);
    double const delta_k = delta_km1 + delta_kp1;
    O_mat(k, k) += coeff / (delta_km1 * delta_km1 * delta_k * delta_k);
    O_mat(k + 2, k + 2) += coeff / (delta_kp1 * delta_kp1 * delta_k * delta_k);
    O_mat(k + 1, k + 1) +=
        coeff / (delta_km1 * delta_km1 * delta_kp1 * delta_kp1);
    O_mat(k, k + 2) += coeff / (delta_km1 * delta_kp1 * delta_k * delta_k);
    O_mat(k + 2, k) += coeff / (delta_km1 * delta_kp1 * delta_k * delta_k);
    O_mat(k, k + 1) += -coeff / (delta_km1 * delta_km1 * delta_kp1 * delta_k);
    O_mat(k + 1, k) += -coeff / (delta_km1 * delta_km1 * delta_kp1 * delta_k);
    O_mat(k + 1, k + 2) +=
        -coeff / (delta_km1 * delta_kp1 * delta_kp1 * delta_k);
    O_mat(k + 2, k + 1) +=
        -coeff / (delta_km1 * delta_kp1 * delta_kp1 * delta_k);
  }

#ifdef EXT_DEBUG
  using triqs::utility::is_zero;
  auto O_mat_view = O_mat(r_K, r_K);
  if(!is_zero(max_element(abs(O_mat_view - transpose(O_mat_view))), 1e-10))
    TRIQS_RUNTIME_ERROR << "Matrix O is not symmetric, O = " << O_mat_view;
#endif
}

template <typename KernelType>
void update_consistent_constraints<KernelType>::optimize_heights() {
  auto K = long(proposed_conf.size());

  // Compute distances between centers of adjacent rectangles
  for(auto k : r_K1) {
    delta_c(k) = proposed_conf[k + 1].center - proposed_conf[k].center;
  }

  Q0(r_K) = 0;
  Q1(r_K1) = Q1_init;
  Q2(r_K2) = Q2_init;

  auto QF_mat_view =
      matrix_view<double, nda::F_layout>({K, K}, QF_mat_data.data());

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
    QF_mat_view() = chi2_mat(r_K, r_K) + O_mat(r_K, r_K);

#ifdef EXT_DEBUG
    double QF = NAN;
#endif
    try {
#ifdef EXT_DEBUG
      QF =
#endif
          worker(QF_mat_view,
                 chi2_rhs(r_K),
                 widths_matrix_view,
                 norm_vector_view,
                 heights(r_K));
    } catch(triqs::runtime_error&) {
#ifdef EXT_DEBUG
      std::cerr << "Cholesky factorization has failed\n";
#endif
      break;
    }

    double norm_diff_max = sum(abs((make_array_view(heights(r_K)) -
                                    make_array_view(heights_prev(r_K))) *
                                   widths(r_K))) /
                           norm;

#ifdef EXT_DEBUG
    double O = nda::blas::dot(heights(r_K), O_mat(r_K, r_K) * heights(r_K));
    double chi2 = QF + chi2_const - O;

    std::cerr << "  Iteration " << (iter + 1) << " out of " << max_iter << ": "
              << "sum(|Δnorm_k|) / norm = " << norm_diff_max
              << ", χ^2 = " << chi2 << ", O/χ^2 = " << (O / chi2) << '\n';
#endif

    if(norm_diff_max < rect_norm_variation_tol) {
#ifdef EXT_DEBUG
      nda::array<double, 1> rectangle_norms =
          make_array_view(heights(r_K)) * widths(r_K);
      std::cerr << "Convergence reached, rectangle norms = " << rectangle_norms
                << '\n';
#endif
      break;
    }

    // Update amplitude regularization parameters
    for(auto k : r_K) Q0(k) = (heights(k) < 0) ? Q0_max : (Q0(k) / Q0_divisor);

    // Update 1st derivative regularization parameters
    finite_diff_forward_dx(heights(r_K), delta_c(r_K1), conf_1st_der(r_K1));

    double Q1_limit = Q12_limiter * min_element(Q1(r_K1));
    for(auto k : r_K1) {
      double d = std::abs(conf_1st_der(k));
      if(d > C / Q1(k))
        Q1(k) = C / d;
      else
        Q1(k) *= Q12_increase_coeff;
      if(Q1(k) > Q1_limit) Q1(k) = Q1_limit;
    }

    // Update 2nd derivative regularization parameters
    finite_diff_2_symm_dx(heights(r_K), delta_c(r_K1), conf_2nd_der(r_K2));

    double Q2_limit = Q12_limiter * min_element(Q2(r_K2));
    for(auto k : r_K2) {
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
    std::cerr << "Maximum number of iterations reached, heights = "
              << heights(r_K) << '\n';
  }
#endif
}

template <typename KernelType>
double update_consistent_constraints<KernelType>::get_acceptance_rate() const {
  return double(n_accepted) / double(n_proposed);
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(update_consistent_constraints)

} // namespace som
