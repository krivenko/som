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

#include <boost/math/special_functions/digamma.hpp>

#include <triqs/utility/numeric_ops.hpp>

#include "fermiongf_imtime.hpp"

namespace som {

using namespace triqs::gfs;

//////////////////////////////////////////
// kernel<FermionGf, imtime>::evaluator //
//////////////////////////////////////////

kernel<FermionGf, imtime>::evaluator::evaluator(
    mesh_type const& mesh,
    mesh_type::mesh_point_t const& tau) {
  using boost::math::digamma;
  using std::atan;
  using std::exp;
  using std::expm1;
  using std::log;
  using std::log1p;

  // Limit for spline interpolation
  static const double x0 = -2 * std::log(tolerance);

  size_t i = tau.linear_index();
  size_t s = mesh.size();
  double beta_ = mesh.domain().beta;
  alpha = double(tau) / beta_;
  double dx = x0 / (n_spline_knots - 1);

  vector<double> spline_m_knots(n_spline_knots), spline_p_knots(n_spline_knots);

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

  if(i < s / 2) { // \alpha < 1/2
    if(i == 0) {  // \alpha = 0
      alpha_case = zero;
      nda::for_each(spline_m_knots.shape(), [&spline_m_knots, dx](int xi) {
        spline_m_knots(xi) = log1p(exp(-x0 + dx * xi)) - log(2.0);
      });
      nda::for_each(spline_p_knots.shape(), [&spline_p_knots, dx](int xi) {
        spline_p_knots(xi) = log1p(exp(-dx * xi)) - log(2.0);
      });
      tail_coeff = -1 / beta_;
    } else { // \alpha \in (0;1/2)
      alpha_case = small;
      // Fill spline_m_knots
      double shift =
          0.5 * (digamma(1 - 0.5 * alpha) - digamma(0.5 - 0.5 * alpha));
      for(auto xi : range(0, n_spline_knots - 1)) {
        double x = -x0 + dx * double(xi);
        spline_m_knots[xi] =
            aux_sum([this](int n) { return n + 1 - alpha; }, x) - shift;
      }
      spline_m_knots[n_spline_knots - 1] = 0;
      // Fill spline_p_knots
      shift = 0.5 * (digamma(1 + 0.5 * alpha) - digamma(0.5 + 0.5 * alpha));
      spline_p_knots[0] = 0;
      for(auto xi : range(1, n_spline_knots)) {
        double x = dx * double(xi);
        spline_p_knots[xi] =
            -aux_sum([this](int n) { return -(n + 1 + alpha); }, x) - shift;
      }
      tail_coeff = 1 / (beta_ * alpha);
    }
  } else if(i >= s - s / 2) { // \alpha > 1/2
    if(i == s - 1) {          // \alpha = 1
      alpha_case = one;
      nda::for_each(spline_m_knots.shape(), [&spline_m_knots, dx](int xi) {
        spline_m_knots(xi) = -log1p(exp(-x0 + dx * xi)) + log(2.0);
      });
      nda::for_each(spline_p_knots.shape(), [&spline_p_knots, dx](int xi) {
        spline_p_knots(xi) = -log1p(exp(-dx * xi)) + log(2.0);
      });
      tail_coeff = -1 / beta_;
    } else { // \alpha \in (1/2;1)
      alpha_case = big;
      // Fill spline_m_knots
      double shift =
          0.5 * (digamma(1.5 - 0.5 * alpha) - digamma(1 - 0.5 * alpha));
      for(auto xi : range(0, n_spline_knots - 1)) {
        double x = -x0 + dx * double(xi);
        spline_m_knots[xi] =
            -aux_sum([this](int n) { return n + 2 - alpha; }, x) + shift;
      }
      spline_m_knots[n_spline_knots - 1] = 0;
      // Fill spline_p_knots
      shift = 0.5 * (digamma(0.5 + 0.5 * alpha) - digamma(0.5 * alpha));
      spline_p_knots[0] = 0;
      for(auto xi : range(1, n_spline_knots)) {
        double x = dx * double(xi);
        spline_p_knots[xi] =
            aux_sum([this](int n) { return -(n + alpha); }, x) + shift;
      }
      tail_coeff = -1 / (beta_ * (1 - alpha));
    }
  } else { // \alpha = 1/2
    alpha_case = half;
    nda::for_each(spline_m_knots.shape(), [&spline_m_knots, dx](int xi) {
      spline_m_knots(xi) = 2 * atan(exp(0.5 * (-x0 + dx * xi))) - M_PI / 2;
    });
    nda::for_each(spline_p_knots.shape(), [&spline_p_knots, dx](int xi) {
      spline_p_knots(xi) = -2 * atan(exp(-0.5 * (dx * xi))) + M_PI / 2;
    });
  }

  spline_m_knots /= -beta_;
  spline_p_knots /= -beta_;
  spline_m = regular_spline(-x0, .0, spline_m_knots);
  spline_p = regular_spline(.0, x0, spline_p_knots);
}

double kernel<FermionGf, imtime>::evaluator::operator()(double x) const {
  using std::expm1;
  switch(alpha_case) {
    case zero: return (x > 0 ? spline_p(x) + tail_coeff * x : spline_m(x));
    case small:
      return (x > 0 ? spline_p(x) + tail_coeff * expm1(-alpha * x)
                    : spline_m(x));
    case half: return (x > 0 ? spline_p(x) : spline_m(x));
    case big:
      return (x >= 0 ? spline_p(x)
                     : spline_m(x) + tail_coeff * expm1((1 - alpha) * x));
    case one: return (x >= 0 ? spline_p(x) : spline_m(x) + tail_coeff * x);
    default: TRIQS_RUNTIME_ERROR << "Internal error: invalid alpha_case";
  }
}

///////////////////////////////
// kernel<FermionGf, imtime> //
///////////////////////////////

kernel<FermionGf, imtime>::kernel(mesh_type const& mesh)
   : kernel_base(mesh.size()), beta(mesh.domain().beta), mesh(mesh) {

  lambdas.reserve(mesh.size());
  for(auto tau : mesh) lambdas.emplace_back(mesh, tau);
}

// Apply to a rectangle
void kernel<FermionGf, imtime>::apply(rectangle const& rect,
                                      result_view_type res) const {

  double x1 = beta * rect.left();
  double x2 = beta * rect.right();

  for(auto tau : mesh) {
    auto const& lambda = lambdas[tau.linear_index()];
    res(tau.linear_index()) = rect.height * (lambda(x2) - lambda(x1));
  }
}

std::ostream& operator<<(std::ostream& os,
                         kernel<FermionGf, imtime> const& kern) {
  os << R"(A(ϵ) -> G(τ), )";
  os << R"(β = )" << kern.beta << ", " << kern.mesh.size() << R"( τ-points)";
  return os;
}

} // namespace som
