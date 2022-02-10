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

#include <boost/math/special_functions/bessel.hpp>

#include <triqs/utility/numeric_ops.hpp>

#include "../numerics/simpson.hpp"
#include "zerotemp_legendre.hpp"

namespace som {

using namespace triqs::gfs;

// Make coefficient of a Bessel polynomial a_k(l+1/2)
static double make_a(int k, int l) {
  double a = 1;
  for(int i = 1; i <= k; ++i) {
    double t = l - k + 2 * i;
    a *= (t - 1) * t / double(2 * i);
  }
  return a;
}

///////////////////////////////////////////
// kernel<ZeroTemp, legendre>::evaluator //
///////////////////////////////////////////

kernel<ZeroTemp, legendre>::evaluator::evaluator(int l, double x0_start)
   : sqrt_pref((l % 2 ? 1 : -1) * std::sqrt(2 * l + 1)) {

  // Integrand, 2 i_l(x) exp(-x)
  auto integrand = [l](double x) {
    if(x == 0) return (l == 0 ? 2.0 : 0.0);
    double val = boost::math::cyl_bessel_i(l + 0.5, x);
    return val * std::sqrt(M_PI / (2 * x)) * 2 * std::exp(-x);
  };

  vector<double> tail_coeffs(l + 2);
  tail_coeffs[0] = 0;
  for(int k = 0; k <= l; ++k)
    tail_coeffs[k + 1] = ((k % 2) ? -1 : 1) * make_a(k, l);
  polynomial<> integrand_tail(tail_coeffs);

  // Search for the low-energy/high-energy boundary
  x0 = x0_start;
  using triqs::utility::is_zero;
  while(!is_zero(integrand(x0) - integrand_tail(1 / x0), tolerance)) {
    x0 += x0_step;
    // http://en.cppreference.com/w/cpp/numeric/math/exp
    static const int max_stable_exp_arg = 708;
    if(x0 >= max_stable_exp_arg)
      TRIQS_RUNTIME_ERROR << "kernel<ZeroTemp,legendre>: l = " << l
                          << " is too large and may cause numerical overflows.";
  }

  // Fill high_energy_pol
  vector<double> int_tail_coeffs(l + 1);
  int_tail_coeffs[0] = 0;
  for(int k = 1; k <= l; ++k) int_tail_coeffs[k] = -tail_coeffs[k + 1] / k;
  high_energy_pol = polynomial<>(int_tail_coeffs);

  // Fill low_energy_spline
  auto spline_knots = primitive(integrand, .0, x0, n_spline_knots, tolerance);
  low_energy_spline = regular_spline(0, x0, spline_knots);

  // Set low_energy_x0 and high_energy_x0
  low_energy_x0 = low_energy_spline(x0);
  high_energy_x0 = std::log(x0) + high_energy_pol(1 / x0);
}

double kernel<ZeroTemp, legendre>::evaluator::operator()(double x) const {
  double val = (x <= x0) ? low_energy_spline(x)
                         : (low_energy_x0 - high_energy_x0 +
                            (std::log(x) + high_energy_pol(1 / x)));
  return sqrt_pref * val;
}

////////////////////////////////
// kernel<ZeroTemp, legendre> //
////////////////////////////////

kernel<ZeroTemp, legendre>::kernel(mesh_type const& mesh)
   : kernel_base(mesh.size()), mesh(mesh), tau_max(mesh.domain().beta) {
  evaluators.reserve(mesh.size());

  double x0 = x0_start_l0;
  for(auto l : mesh) {
    evaluators.emplace_back(l, x0);
    x0 = evaluators.back().x0;
  }
}

void kernel<ZeroTemp, legendre>::apply(rectangle const& rect,
                                       result_type& res) const {

  double e1 = rect.left();
  double e2 = rect.right();

  for(auto l : mesh)
    res(l.linear_index()) = rect.height * (Lambda(l, e2) - Lambda(l, e1));
}

double kernel<ZeroTemp, legendre>::Lambda(int l, double Omega) const {
  return evaluators[l](Omega * tau_max / 2);
}

std::ostream& operator<<(std::ostream& os,
                         kernel<ZeroTemp, legendre> const& kern) {
  os << R"(A(\epsilon) -> G_{T=0}(l), )";
  os << R"(\tau_{max} = )" << kern.tau_max << ", " << kern.mesh.size()
     << " Legendre coefficients.";
  return os;
}

} // namespace som
