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

#include <iostream>
#include <vector>

#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>

#include "../numerics/polynomial.hpp"
#include "../numerics/spline.hpp"

#include "base.hpp"

namespace som {

// Kernel: fermionic GF, Legendre basis
template <>
class kernel<BosonCorr, triqs::gfs::legendre>
   : public kernel_base<kernel<BosonCorr, triqs::gfs::legendre>,
                        triqs::arrays::array<double, 1>> {

  // Tolerance levels for function evaluation
  static constexpr double tolerance = 1e-14;
  // Number of energy points for spline interpolation
  static constexpr int n_spline_knots = 15001;
  // Starting low-energy/high-energy boundary value for l=0
  static constexpr double x0_start_l0 = 15.0;
  // Step used in x0 boundary search
  static constexpr double x0_step = 1.0;

  // Evaluator object
  // operator(z) returns
  // (2/(\pi*\beta)) * sqrt(2*l+1) * \int_0^z i_l(x) / sinh(x) dx
  struct evaluator {
    double x0; // Boundary between low-energy and high-energy regions
    regular_spline low_energy_spline; // Spline approximation on [0;x0]
    polynomial<> high_energy_pol;     // Polynomial approximation on ]x0;\infty[
    double log_coeff;                 // -a_1(l+1/2) = -l(l+1)/2
    double low_energy_x0;             // low_energy_spline(x0)
    double high_energy_x0; // x - log_coeff * log(x0) + high_energy_pol(1/x0)
    double pref;           // (2/(\pi*\beta)) * sqrt(2*l+1) prefactor

    evaluator(int l, double x0_start, double beta);

    double operator()(double x) const;
  };

  // Evaluator objects, one object per l
  std::vector<evaluator> evaluators;

  // Integrated kernel \Lambda(l,\Omega)
  double Lambda(int l, double Omega) const;

public:
  using result_type = triqs::arrays::array<double, 1>;
  using mesh_type = triqs::gfs::gf_mesh<triqs::gfs::legendre>;
  constexpr static observable_kind kind = BosonCorr;

  const double beta;    // Inverse temperature
  const mesh_type mesh; // Legendre coefficients mesh

  explicit kernel(mesh_type const& mesh);

  // Apply to a rectangle
  void apply(rectangle const& rect, result_type& res) const;

  friend std::ostream& operator<<(std::ostream& os, kernel const& kern);
};

} // namespace som
