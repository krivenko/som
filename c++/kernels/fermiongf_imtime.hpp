/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <vector>
#include <cmath>
#include <triqs/gfs.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include "base.hpp"
#include "../numerics/spline.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

// Kernel: fermionic GF, imaginary time
template<> class kernel<FermionGf,imtime> :
           public kernel_base<kernel<FermionGf,imtime>, array<double,1>> {

public:

 using result_type = array<double,1>;
 using mesh_type = gf_mesh<imtime>;
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
  enum {zero, small, half, big, one} alpha_case;

  // Coefficient of the tail function
  double tail_coeff;

  // S(x) = \sum_{n=0}^{+\infty} \frac{(-1)^n \exp(d(n) x)}{d(n)}
  template<typename F> static double aux_sum(F && d, double x) {
   using triqs::utility::is_zero;
   double val = 0;
   for(int n = 0;;++n) {
    double dd = d(n);
    double t = (n % 2 ? -1 : 1) * std::exp(dd*x) / dd;
    if(is_zero(t, tolerance)) break;
    val += t;
   }
   return val;
  }

  evaluator(mesh_type const& mesh, mesh_type::mesh_point_t const& tau) {
   using std::exp;
   using std::expm1;
   using std::log;
   using std::log1p;
   using std::atan;
   using boost::math::digamma;

   // Limit for spline interpolation
   static const double x0 = -2*std::log(tolerance);

   int i = tau.index();
   int s = mesh.size();
   double beta = mesh.domain().beta;
   alpha = double(tau) / beta;
   double dx = x0 / (n_spline_knots - 1);

   vector<double> spline_m_knots(n_spline_knots), spline_p_knots(n_spline_knots);

   if(i < s/2) {             // \alpha < 1/2
    if(i == 0) {             // \alpha = 0
     alpha_case = zero;
     assign_foreach(spline_m_knots, [dx](int xi){ return log1p(exp(-x0 + dx*xi)) - log(2.0); });
     assign_foreach(spline_p_knots, [dx](int xi){ return log1p(exp(-dx*xi)) - log(2.0); });
     tail_coeff = -1/beta;
    } else {                 // \alpha \in (0;1/2)
     alpha_case = small;
     // Fill spline_m_knots
     double shift = 0.5*(digamma(1-0.5*alpha) - digamma(0.5-0.5*alpha));
     for(int xi : range(0, n_spline_knots-1)) {
      double x = -x0+dx*xi;
      spline_m_knots[xi] = aux_sum([this](int n){return n+1-alpha; }, x) - shift;
     }
     spline_m_knots[n_spline_knots-1] = 0;
     // Fill spline_p_knots
     shift = 0.5*(digamma(1+0.5*alpha) - digamma(0.5+0.5*alpha));
     spline_p_knots[0] = 0;
     for(int xi : range(1, n_spline_knots)) {
      double x = dx*xi;
      spline_p_knots[xi] = -aux_sum([this](int n){return -(n+1+alpha); }, x) - shift;
     }
     tail_coeff = 1/(beta * alpha);
    }
   } else if(i >= s - s/2) { // \alpha > 1/2
    if(i == s-1) {           // \alpha = 1
     alpha_case = one;
     assign_foreach(spline_m_knots, [dx](int xi){ return -log1p(exp(-x0 + dx*xi)) + log(2.0); });
     assign_foreach(spline_p_knots, [dx](int xi){ return -log1p(exp(-dx*xi)) + log(2.0); });
     tail_coeff = -1/beta;
    } else {                 // \alpha \in (1/2;1)
     alpha_case = big;
     // Fill spline_m_knots
     double shift = 0.5*(digamma(1.5-0.5*alpha) - digamma(1-0.5*alpha));
     for(int xi : range(0, n_spline_knots-1)) {
      double x = -x0+dx*xi;
      spline_m_knots[xi] = -aux_sum([this](int n){return n+2-alpha; }, x) + shift;
     }
     spline_m_knots[n_spline_knots-1] = 0;
     // Fill spline_p_knots
     shift = 0.5*(digamma(0.5+0.5*alpha) - digamma(0.5*alpha));
     spline_p_knots[0] = 0;
     for(int xi : range(1, n_spline_knots)) {
      double x = dx*xi;
      spline_p_knots[xi] = aux_sum([this](int n){return -(n+alpha); }, x) + shift;
     }
     tail_coeff = -1/(beta * (1-alpha));
    }
   } else {                  // \alpha = 1/2
    alpha_case = half;
    assign_foreach(spline_m_knots, [dx](int xi){ return 2*atan(exp(0.5*(-x0+dx*xi))) - M_PI/2; });
    assign_foreach(spline_p_knots, [dx](int xi){ return -2*atan(exp(-0.5*(dx*xi))) + M_PI/2; });
   }

   spline_m_knots /= -beta; spline_p_knots /= -beta;
   spline_m = regular_spline(-x0, .0, spline_m_knots);
   spline_p = regular_spline( .0, x0, spline_p_knots);
  }

  double operator()(double x) const {
   using std::expm1;
   switch(alpha_case) {
    case zero: return (x > 0 ? spline_p(x) + tail_coeff * x : spline_m(x));
    case small: return (x > 0 ? spline_p(x) + tail_coeff * expm1(-alpha*x) : spline_m(x));
    case half: return (x > 0 ? spline_p(x) : spline_m(x));
    case big: return (x >= 0 ? spline_p(x) : spline_m(x) + tail_coeff * expm1((1-alpha)*x));
    case one: return (x >= 0 ? spline_p(x) : spline_m(x) + tail_coeff * x);
    default: TRIQS_RUNTIME_ERROR << "Internal error: invalid alpha_case";
   }
  }

 };

 // List of integrated kernels for all \alpha
 std::vector<evaluator> lambdas;

public:

 const double beta;          // Inverse temperature
 const mesh_type mesh;       // Matsubara time mesh

 kernel(mesh_type const& mesh) :
  kernel_base(mesh.size()), mesh(mesh), beta(mesh.domain().beta) {

  lambdas.reserve(mesh.size());
  for(auto tau : mesh) lambdas.emplace_back(mesh, tau);
 }

 // Apply to a rectangle
 void apply(rectangle const& rect, result_type & res) const {

  double x1 = beta*(rect.center - rect.width/2);
  double x2 = beta*(rect.center + rect.width/2);

  for(auto tau : mesh) {
   auto const& lambda = lambdas[tau.index()];
   res(tau.index()) = rect.height * (lambda(x2) - lambda(x1));
  }
 }

 friend std::ostream & operator<<(std::ostream & os, kernel const& kern) {
  os << "A(\\epsilon) -> G(\\tau), ";
  os << "\\beta = " << kern.beta << ", " << kern.mesh.size() << " \\tau-points.";
  return os;
 }

};

}
