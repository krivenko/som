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
#include <cmath>

#include <boost/math/special_functions/digamma.hpp>

#include <triqs/utility/numeric_ops.hpp>

#include "fermiongfsymm_imtime.hpp"

namespace som {

using namespace triqs::gfs;

//////////////////////////////////////////////
// kernel<FermionGfSymm, imtime>::evaluator //
//////////////////////////////////////////////

kernel<FermionGfSymm, imtime>::evaluator::evaluator(
    mesh_type const& mesh,
    mesh_type::mesh_point_t const& tau) {
  using boost::math::digamma;

  // Limit for spline interpolation
  static const double x0 = -2 * std::log(tolerance);

  size_t i = tau.data_index();
  size_t s = mesh.size();
  double beta_ = mesh.beta();
  alpha = double(tau) / beta_;
  double dx = x0 / (n_spline_knots - 1);

  vector<double> spline_knots(n_spline_knots);

  if(i == 0 || i == s - 1) { // \alpha = 0, 1
    alpha_case = edge;
    spline_knots() = 0;
    tail_coeff1 = -0.5 / beta_;
  } else if(s % 2 == 1 && i == s / 2) { // \alpha = 1/2
    alpha_case = half;
    nda::for_each(spline_knots.shape(), [&spline_knots, dx](int xi) {
      if(xi == 0) {
        spline_knots(xi) = 0;
        return;
      }
      spline_knots(xi) = 4 * std::atan(std::tanh(dx * xi / 4));
    });
  } else { // \alpha \in (0;1/2)\cup(1/2;1)
    // S(x) = \sum_{n=0}^{+\infty} \frac{(-1)^n \exp(d(n) x)}{d(n)}
    auto aux_sum = [](auto d, double x) -> double {
      using triqs::utility::is_zero;
      double val = 0;
      for(int n = 0;; ++n) {
        double dd = d(n);
        double t = ((n % 2) ? -1 : 1) * std::exp(dd * x) / dd;
        if(is_zero(t, tolerance)) break;
        val += t;
      }
      return val;
    };
    alpha_case = other;
    spline_knots[0] = 0;
    double shift =
        0.5 * (digamma(1 + 0.5 * alpha) - digamma(0.5 + 0.5 * alpha)) +
        0.5 * (digamma(1.5 - 0.5 * alpha) - digamma(1.0 - 0.5 * alpha));
    for(auto xi : range(1, n_spline_knots)) {
      double x = dx * double(xi);
      spline_knots[xi] =
          -aux_sum([this](int n) { return -(n + 1 + alpha); }, x) -
          aux_sum([this](int n) { return -(n + 2 - alpha); }, x) - shift;
    }

    tail_coeff1 = 1 / (2 * alpha * beta_);
    tail_coeff2 = 1 / (2 * (1 - alpha) * beta_);
  }

  spline_knots /= -(2 * beta_);
  spline_ = regular_spline(.0, x0, spline_knots);
}

double kernel<FermionGfSymm, imtime>::evaluator::operator()(double x) const {
  if(x < 0) return -operator()(-x);

  using std::expm1;
  switch(alpha_case) {
    case edge: return tail_coeff1 * x;
    case other:
      return spline_(x) + tail_coeff1 * expm1(-alpha * x) +
             tail_coeff2 * expm1(-(1 - alpha) * x);
    case half: return spline_(x);
    default: TRIQS_RUNTIME_ERROR << "Internal error: invalid alpha_case";
  }
}

///////////////////////////////////
// kernel<FermionGfSymm, imtime> //
///////////////////////////////////

kernel<FermionGfSymm, imtime>::kernel(mesh_type const& mesh)
   : kernel_base(mesh.size()), beta(mesh.beta()), mesh(mesh) {

  lambdas.reserve(mesh.size());
  for(auto tau : mesh) lambdas.emplace_back(mesh, tau);
}

// Apply to a rectangle
void kernel<FermionGfSymm, imtime>::apply(rectangle const& rect,
                                          result_view_type res) const {

  double x1 = beta * rect.left();
  double x2 = beta * rect.right();

  for(auto tau : mesh) {
    auto const& lambda = lambdas[tau.data_index()];
    res(tau.data_index()) = rect.height * (lambda(x2) - lambda(x1));
  }
}

std::ostream& operator<<(std::ostream& os,
                         kernel<FermionGfSymm, imtime> const& kern) {
  os << R"(A(ϵ) -> G_{sym}(τ), )";
  os << R"(β = )" << kern.beta << ", " << kern.mesh.size() << R"( τ-points)";
  return os;
}

} // namespace som
