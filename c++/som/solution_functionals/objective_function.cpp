/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <algorithm>

#include <nda/linalg.hpp>

#include <triqs/utility/exceptions.hpp>

#include "objective_function.hpp"

namespace som {

template <typename KernelType>
objective_function<KernelType>::objective_function(
    KernelType const& kern,
    rhs_type const& rhs,
    error_bars_type const& error_bars)
   : kern(kern), rhs(rhs), sigma2(abs2(error_bars)), tmp(rhs.size()) {}

template <typename KernelType>
objective_function<KernelType>::objective_function(
    KernelType const& kern,
    rhs_type const& rhs,
    cov_matrix_type const& cov_matrix,
    double filtering_level)
   : kern(kern), rhs(rhs), tmp(rhs.size()) {
  auto [ev, vecs] = nda::linalg::eigenelements(cov_matrix);

  auto first_positive_sigma2 = std::distance(
      std::begin(ev), std::upper_bound(std::begin(ev), std::end(ev), 0));
  if(first_positive_sigma2 == ev.size())
    TRIQS_RUNTIME_ERROR << "objective_function: Covariance matrix does not "
                           "have a single positive eigenvalue";

  auto r = nda::range(first_positive_sigma2, ev.size());
  sigma2 = ev(r);
  sigma2() += filtering_level * filtering_level;
  U_dagger = dagger(vecs(nda::range::all, r));
}

template <typename KernelType>
inline double objective_function<KernelType>::call_impl() const {
  if(U_dagger) {
    auto r = nda::range(sigma2.size());
    tmp(r) = (*U_dagger) * nda::vector_view<rhs_scalar_type>(tmp);
    return sum(abs2(tmp(r)) / sigma2) / double(sigma2.size());
  } else {
    return sum(abs2(tmp) / sigma2) / double(rhs.size());
  }
}

template <typename KernelType>
double
objective_function<KernelType>::operator()(configuration const& c) const {
  tmp() = rhs - kern(c);
  return call_impl();
}

template <typename KernelType>
double
objective_function<KernelType>::operator()(config_update const& cu) const {
  tmp() = rhs - kern(cu);
  return call_impl();
}

INSTANTIATE_CLASS_FOR_EACH_KERNEL(objective_function)

} // namespace som
