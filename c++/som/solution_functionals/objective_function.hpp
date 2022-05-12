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

#include <optional>

#include <nda/nda.hpp>

#include <som/config_update.hpp>
#include <som/configuration.hpp>
#include <som/kernels/all.hpp>

namespace som {

// Objective function \chi^2
template <typename KernelType> class objective_function {

public:
  using rhs_type = typename KernelType::result_type;
  using mesh_type = typename KernelType::mesh_type;

  using rhs_scalar_type = typename rhs_type::value_type;

  using error_bars_type = rhs_type;
  using cov_matrix_type = nda::matrix<rhs_scalar_type>;

private:
  // Integral kernel
  KernelType const& kern;
  // The right-hand side of the Fredholm integral equation
  rhs_type const& rhs;
  // Shifted positive eigenvalues of the covariance matrix,
  // \sigma_n^2 + (filtration_level)^2.
  nda::array<double, 1> sigma2;

  // Matrix U^\dagger from the eigendecomposition of the covariance matrix,
  // cov = U diag(\sigma_n^2) U^\dagger.
  // The number of rows of U^\dagger equals the number of retained eigenvalues
  // (size of sigma2).
  std::optional<cov_matrix_type> U_dagger;

  // Stores (rhs - kern(c)) and U^\dagger*(rhs - kern(c))
  mutable rhs_type tmp;

  // Common part of both operator() overloads.
  double call_impl() const;

public:
  objective_function(KernelType const& kern,
                     rhs_type const& rhs,
                     error_bars_type const& error_bars);

  objective_function(KernelType const& kern,
                     rhs_type const& rhs,
                     cov_matrix_type const& cov_matrix,
                     double filtration_level = 0);

  double operator()(configuration const& c) const;
  double operator()(config_update const& cu) const;

  [[nodiscard]] KernelType const& get_kernel() const { return kern; }
  [[nodiscard]] rhs_type const& get_rhs() const { return rhs; }
  [[nodiscard]] nda::array<double, 1> const& get_sigma2() const {
    return sigma2;
  }
  [[nodiscard]] std::optional<cov_matrix_type> get_U_dagger() const {
    return U_dagger;
  }
};

EXTERN_TEMPLATE_CLASS_FOR_EACH_KERNEL(objective_function)

} // namespace som
