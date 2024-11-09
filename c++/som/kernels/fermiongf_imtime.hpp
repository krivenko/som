/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko
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

#include <nda/nda.hpp>
#include <triqs/mesh.hpp>

#include "../numerics/spline.hpp"

#include "base.hpp"

namespace som {

// Kernel: fermionic GF, imaginary time
template <>
class kernel<FermionGf, triqs::mesh::imtime>
   : public kernel_base<kernel<FermionGf, triqs::mesh::imtime>,
                        nda::array<double, 1>> {

public:
  using result_type = nda::array<double, 1>;
  using result_view_type = nda::array_view<double, 1>;
  using mesh_type = triqs::mesh::imtime;
  constexpr static observable_kind kind = FermionGf;

private:
  // Integrated kernel \Lambda(\alpha,x)
  struct evaluator {
    // Tolerance levels for function evaluation
    static constexpr double tolerance = 1e-14;
    // Number of x points for Lambda interpolation
    static constexpr int n_spline_knots = 10001;

    // Spline interpolations S^-(x) and S^+(x)
    regular_spline spline_m, spline_p;

    // \alpha = \tau / \beta
    double alpha;

    // Value of \alpha
    // zero: \alpha = 0
    // small: \alpha \in (0;1/2)
    // half: \alpha = 1/2
    // big: \alpha \in (1/2;1)
    // one: \alpha = 1
    enum { zero, small, half, big, one } alpha_case;

    // Coefficient of the tail function
    double tail_coeff;

    evaluator(mesh_type const& mesh, mesh_type::mesh_point_t const& tau);

    double operator()(double x) const;
  };

  // List of integrated kernels for all \alpha
  std::vector<evaluator> lambdas;

  const double beta;    // Inverse temperature
  const mesh_type mesh; // Matsubara time mesh

public:
  explicit kernel(mesh_type const& mesh);

  // Apply to a rectangle
  void apply(rectangle const& rect, result_view_type res) const;

  friend std::ostream& operator<<(std::ostream& os, kernel const& kern);
};

} // namespace som
