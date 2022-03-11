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

#include <boost/math/special_functions/trigamma.hpp>

#include <triqs/utility/numeric_ops.hpp>

#include "../numerics/dilog.hpp"
#include "bosonautocorr_imtime.hpp"

namespace som {

using namespace triqs::gfs;

//////////////////////////////////////////////
// kernel<BosonAutoCorr, imtime>::evaluator //
//////////////////////////////////////////////

kernel<BosonAutoCorr, imtime>::evaluator::evaluator(
    mesh_type const& mesh, mesh_type::mesh_point_t const& tau) {
  using boost::math::trigamma;
  using som::dilog;
  using std::exp;
  using std::log;
  using std::log1p;
  using std::real;
  using std::tanh;

  // Limit for spline interpolation
  static const double x0 = -1.1 * std::log(tolerance);

  size_t i = tau.linear_index();
  size_t s = mesh.size();
  double beta_ = mesh.domain().beta;
  alpha = double(tau) / beta_;
  double dx = x0 / (n_spline_knots - 1);

  vector<double> spline_knots(n_spline_knots);

  if(i == 0 || i == s - 1) { // \alpha = 0, 1
    alpha_case = edge;
    nda::for_each(spline_knots.shape(), [&spline_knots, dx](int xi) {
      if(xi == 0) { spline_knots(xi) = 0; return; }
      double x = dx * xi;
      double expx = exp(-x);
      spline_knots(xi) = -1 + M_PI * M_PI / 3 + 2 * x * log1p(-expx) -
                         2 * real(dilog(expx)) + expx * (1 + x);
    });
    tail_coeff1 = 1 / (M_PI * beta_ * beta_);
  } else if(s % 2 == 1 && i == s / 2) { // \alpha = 1/2
    alpha_case = half;
    nda::for_each(spline_knots.shape(), [&spline_knots, dx](int xi) {
      if(xi == 0) { spline_knots(xi) = 0; return; }
      double x = dx * xi;
      double expx2 = exp(-x / 2);
      spline_knots(xi) = -8 + M_PI * M_PI + 4 * expx2 * (2 + x) + 2 * x * log(tanh(x / 4)) -
             8 * real(dilog(expx2)) + 2 * real(dilog(expx2 * expx2));
    });
    tail_coeff1 = 4 / (M_PI * beta_ * beta_);
  } else { // \alpha \in (0;1/2)\cup(1/2;1)
    // S(x) = \sum_{n=0}^{+\infty} \frac{\exp(d(n) x)[1-xd(n)]}{d^2(n)}
    auto aux_sum = [](auto d, double x) -> double {
      using triqs::utility::is_zero;
      double val = 0;
      for(int n = 0;; ++n) {
        double dd = d(n);
        double t = std::exp(dd * x) * (1 - x * dd) / (dd * dd);
        if(is_zero(t, tolerance)) break;
        val += t;
      }
      return val;
    };
    alpha_case = other;
    double shift = trigamma(1 + alpha) + trigamma(2 - alpha);
    spline_knots[0] = 0;
    for(int xi : range(1, n_spline_knots)) {
      double x = dx * xi;
      spline_knots[xi] =
          -aux_sum([this](int n) { return -(n + 1 + alpha); }, x) -
          aux_sum([this](int n) { return -(n + 2 - alpha); }, x) + shift;
    }
    tail_coeff1 = 1 / ((M_PI * beta_ * beta_) * alpha * alpha);
    tail_coeff2 = 1 / ((M_PI * beta_ * beta_) * (1 - alpha) * (1 - alpha));
  }

  spline_knots /= (M_PI * beta_ * beta_);
  spline_ = regular_spline(.0, x0, spline_knots);
}

double kernel<BosonAutoCorr, imtime>::evaluator::operator()(double x) const {
  using std::exp;
  switch(alpha_case) {
    case edge:
      return spline_(x) + tail_coeff1 * (1 + x * x / 2 - exp(-x) * (1 + x));
    case half: return spline_(x) + tail_coeff1 * (2 - exp(-x / 2) * (2 + x));
    case other: {
      double y1 = x * alpha, y2 = x * (1 - alpha);
      return spline_(x) + tail_coeff1 * (1 - exp(-y1) * (1 + y1)) +
             tail_coeff2 * (1 - exp(-y2) * (1 + y2));
    }
    default: TRIQS_RUNTIME_ERROR << "Internal error: invalid alpha_case";
  }
}

///////////////////////////////////
// kernel<BosonAutoCorr, imtime> //
///////////////////////////////////

kernel<BosonAutoCorr, imtime>::kernel(mesh_type const& mesh)
   : kernel_base(mesh.size()), beta(mesh.domain().beta), mesh(mesh) {

  lambdas.reserve(mesh.size());
  for(auto tau : mesh) lambdas.emplace_back(mesh, tau);
}

// Apply to a rectangle
void kernel<BosonAutoCorr, imtime>::apply(rectangle const& rect,
                                          result_type& res) const {

  double x1 = beta * rect.left();
  double x2 = beta * rect.right();

  for(auto tau : mesh) {
    auto const& lambda = lambdas[tau.linear_index()];
    res(tau.linear_index()) = rect.height * (lambda(x2) - lambda(x1));
  }
}

std::ostream& operator<<(std::ostream& os,
                         kernel<BosonAutoCorr, imtime> const& kern) {
  os << R"(A(ϵ) -> χ_{sym}(τ), )";
  os << R"(β = )" << kern.beta << ", " << kern.mesh.size() << R"( τ-points.)";
  return os;
}

} // namespace som
