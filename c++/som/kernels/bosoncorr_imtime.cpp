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
#include "bosoncorr_imtime.hpp"

namespace som {

using namespace triqs::gfs;

//////////////////////////////////////////
// kernel<BosonCorr, imtime>::evaluator //
//////////////////////////////////////////

kernel<BosonCorr, imtime>::evaluator::evaluator(
    mesh_type const& mesh, mesh_type::mesh_point_t const& tau) {
  using boost::math::trigamma;
  using som::dilog;
  using std::exp;
  using std::log1p;
  using std::real;
  using std::tanh;

  // Limit for spline interpolation
  static const double x0 = -2.3 * std::log(tolerance);
  static const double pi6 = M_PI * M_PI / 6;

  size_t i = tau.linear_index();
  size_t s = mesh.size();
  double beta_ = mesh.domain().beta;
  alpha = double(tau) / beta_;
  double dx = x0 / (n_spline_knots - 1);

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

  vector<double> spline_m_knots(n_spline_knots), spline_p_knots(n_spline_knots);

  if(i < s / 2) { // \alpha < 1/2
    if(i == 0) {  // \alpha = 0
      alpha_case = zero;
      assign_foreach(spline_m_knots, [dx](int xi) {
        if(xi == n_spline_knots - 1) return .0;
        double x = -x0 + dx * xi;
        double expx = exp(x);
        return real(dilog(expx)) - pi6 + x * log1p(-expx);
      });
      assign_foreach(spline_p_knots, [dx](int xi) {
        if(xi == 0) return .0;
        double x = dx * xi;
        double expx = exp(-x);
        return pi6 - real(dilog(expx)) + x * log1p(-expx);
      });
      tail_coeff = 1 / (M_PI * beta_ * beta_);
    } else { // \alpha \in (0;1/2)
      alpha_case = small;
      // Fill spline_m_knots
      double shift = trigamma(1 - alpha);
      for(int xi : range(0, n_spline_knots - 1)) {
        double x = -x0 + dx * xi;
        spline_m_knots[xi] =
            aux_sum([this](int n) { return n + 1 - alpha; }, x) - shift;
      }
      spline_m_knots[n_spline_knots - 1] = 0;
      // Fill spline_p_knots
      shift = trigamma(1 + alpha);
      spline_p_knots[0] = 0;
      for(int xi : range(1, n_spline_knots)) {
        double x = dx * xi;
        spline_p_knots[xi] =
            -aux_sum([this](int n) { return -(n + 1 + alpha); }, x) + shift;
      }
      tail_coeff = 1 / ((M_PI * beta_ * beta_) * (alpha * alpha));
    }
  } else if(i >= s - s / 2) { // \alpha > 1/2
    if(i == s - 1) {          // \alpha = 1
      alpha_case = one;
      assign_foreach(spline_m_knots, [dx](int xi) {
        if(xi == n_spline_knots - 1) return .0;
        double x = -x0 + dx * xi;
        double expx = exp(x);
        return real(dilog(expx)) - pi6 + x * log1p(-expx);
      });
      assign_foreach(spline_p_knots, [dx](int xi) {
        if(xi == 0) return .0;
        double x = dx * xi;
        double expx = exp(-x);
        return pi6 - real(dilog(expx)) + x * log1p(-expx);
      });
      tail_coeff = 1 / (M_PI * beta_ * beta_);
    } else { // \alpha \in (1/2;1)
      alpha_case = big;
      // Fill spline_m_knots
      double shift = trigamma(2 - alpha);
      for(int xi : range(0, n_spline_knots - 1)) {
        double x = -x0 + dx * xi;
        spline_m_knots[xi] =
            aux_sum([this](int n) { return n + 2 - alpha; }, x) - shift;
      }
      spline_m_knots[n_spline_knots - 1] = 0;
      // Fill spline_p_knots
      shift = trigamma(alpha);
      spline_p_knots[0] = 0;
      for(int xi : range(1, n_spline_knots)) {
        double x = dx * xi;
        spline_p_knots[xi] =
            -aux_sum([this](int n) { return -(n + alpha); }, x) + shift;
      }
      tail_coeff = 1 / ((M_PI * beta_ * beta_) * (1 - alpha) * (1 - alpha));
    }
  } else { // \alpha = 1/2
    alpha_case = half;
    static const double pi2 = M_PI * M_PI / 2;
    assign_foreach(spline_m_knots, [dx](int xi) {
      if(xi == n_spline_knots - 1) return .0;
      double x = -x0 + dx * xi;
      return -pi2 + x * log(tanh(-x / 4)) - real(dilog(exp(x))) +
             4 * real(dilog(exp(x / 2)));
    });
    assign_foreach(spline_p_knots, [dx](int xi) {
      if(xi == 0) return .0;
      double x = dx * xi;
      return pi2 + x * log(tanh(x / 4)) + real(dilog(exp(-x))) -
             4 * real(dilog(exp(-x / 2)));
    });
  }

  spline_m_knots /= (M_PI * beta_ * beta_);
  spline_p_knots /= (M_PI * beta_ * beta_);
  spline_m = regular_spline(-x0, .0, spline_m_knots);
  spline_p = regular_spline(.0, x0, spline_p_knots);
}

double kernel<BosonCorr, imtime>::evaluator::operator()(double x) const {
  using std::exp;
  switch(alpha_case) {
    case zero:
      return (x > 0 ? spline_p(x) + tail_coeff * (x * x / 2) : spline_m(x));
    case small: {
      if(x > 0) {
        double y = x * alpha;
        return spline_p(x) + tail_coeff * (1 - exp(-y) * (1 + y));
      } else
        return spline_m(x);
    }
    case half: return (x > 0 ? spline_p(x) : spline_m(x));
    case big: {
      if(x >= 0)
        return spline_p(x);
      else {
        double y = x * (1 - alpha);
        return spline_m(x) + tail_coeff * (exp(y) * (1 - y) - 1);
      }
    }
    case one:
      return (x >= 0 ? spline_p(x) : spline_m(x) + tail_coeff * (-x * x / 2));
    default: TRIQS_RUNTIME_ERROR << "Internal error: invalid alpha_case";
  }
}

///////////////////////////////
// kernel<BosonCorr, imtime> //
///////////////////////////////

kernel<BosonCorr, imtime>::kernel(mesh_type const& mesh)
   : kernel_base(mesh.size()), beta(mesh.domain().beta), mesh(mesh) {

  lambdas.reserve(mesh.size());
  for(auto tau : mesh) lambdas.emplace_back(mesh, tau);
}

// Apply to a rectangle
void kernel<BosonCorr, imtime>::apply(rectangle const& rect,
                                      result_type& res) const {

  double x1 = beta * rect.left();
  double x2 = beta * rect.right();

  for(auto tau : mesh) {
    auto const& lambda = lambdas[tau.linear_index()];
    res(tau.linear_index()) = rect.height * (lambda(x2) - lambda(x1));
  }
}

std::ostream& operator<<(std::ostream& os,
                         kernel<BosonCorr, imtime> const& kern) {
  os << R"(A(ϵ) -> χ(τ), )";
  os << R"(β = )" << kern.beta << ", " << kern.mesh.size()
     << R"( τ-points.)";
  return os;
}

} // namespace som
