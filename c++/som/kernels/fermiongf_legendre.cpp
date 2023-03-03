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

#include <boost/math/special_functions/bessel.hpp>

#include <triqs/utility/numeric_ops.hpp>

#include "../numerics/simpson.hpp"
#include "fermiongf_legendre.hpp"

namespace som {

using namespace triqs::gfs;

// Make coefficient of a Bessel polynomial a_k(l+1/2)
static double make_a(long k, long l) {
  double a = 1;
  for(long i = 1; i <= k; ++i) {
    auto t = double(l - k + 2 * i);
    a *= (t - 1) * t / double(2 * i);
  }
  return a;
}

////////////////////////////////////////////
// kernel<FermionGf, legendre>::evaluator //
////////////////////////////////////////////

kernel<FermionGf, legendre>::evaluator::evaluator(long l, double x0_start)
   : sqrt_pref(std::sqrt(2 * l + 1)) {

  // Integrand, i_l(x) / cosh(x)
  auto integrand = [l](double x) {
    if(x == 0) return (l == 0 ? 1.0 : 0.0);
    double val = boost::math::cyl_bessel_i(double(l) + 0.5, x);
    return val * std::sqrt(M_PI / (2 * x)) / std::cosh(x);
  };

  vector<double> tail_coeffs(l + 2);
  tail_coeffs[0] = 0;
  for(long k = 0; k <= l; ++k)
    tail_coeffs[k + 1] = ((k % 2) ? -1 : 1) * make_a(k, l);
  polynomial<> integrand_tail(tail_coeffs);

  // Search for the low-energy/high-energy boundary
  x0 = x0_start;
  using triqs::utility::is_zero;
  while(!is_zero(integrand(x0) - integrand_tail(1 / x0), tolerance)) {
    x0 += x0_step;
    // http://en.cppreference.com/w/cpp/numeric/math/cosh
    static const int max_stable_cosh_arg = 710;
    if(x0 >= max_stable_cosh_arg)
      TRIQS_RUNTIME_ERROR << "kernel<FermionGf,legendre>: l = " << l
                          << " is too large and may cause numerical overflows.";
  }

  // Fill high_energy_pol
  vector<double> int_tail_coeffs(l + 1);
  int_tail_coeffs[0] = 0;
  for(long k = 1; k <= l; ++k)
    int_tail_coeffs[k] = -tail_coeffs[k + 1] / double(k);
  high_energy_pol = polynomial<>(int_tail_coeffs);

  // Fill low_energy_spline
  auto spline_knots = primitive(integrand, .0, x0, n_spline_knots, tolerance);
  low_energy_spline = regular_spline(0, x0, spline_knots);

  // Set low_energy_x0 and high_energy_x0
  low_energy_x0 = low_energy_spline(x0);
  high_energy_x0 = std::log(x0) + high_energy_pol(1 / x0);
}

double kernel<FermionGf, legendre>::evaluator::operator()(double x) const {
  double val = (x <= x0) ? low_energy_spline(x)
                         : (low_energy_x0 - high_energy_x0 +
                            (std::log(x) + high_energy_pol(1 / x)));
  return sqrt_pref * val;
}

/////////////////////////////////
// kernel<FermionGf, legendre> //
/////////////////////////////////

kernel<FermionGf, legendre>::kernel(mesh_type const& mesh)
   : kernel_base(mesh.size()), beta(mesh.domain().beta), mesh(mesh) {
  evaluators.reserve(mesh.size());

  double x0 = x0_start_l0;
  for(auto l : mesh) {
    evaluators.emplace_back(l, x0);
    x0 = evaluators.back().x0;
  }
}

void kernel<FermionGf, legendre>::apply(rectangle const& rect,
                                        result_view_type res) const {

  double e1 = rect.left();
  double e2 = rect.right();

  for(auto l : mesh)
    res(l.linear_index()) = rect.height * (Lambda(l, e2) - Lambda(l, e1));
}

double kernel<FermionGf, legendre>::Lambda(long l, double Omega) const {
  return ((l % 2) ? 1.0 : std::copysign(1.0, -Omega)) *
         evaluators[l](std::abs(Omega) * beta / 2);
}

std::ostream& operator<<(std::ostream& os,
                         kernel<FermionGf, legendre> const& kern) {
  os << R"(A(ϵ) -> G(ℓ), )";
  os << R"(β = )" << kern.beta << ", " << kern.mesh.size()
     << " Legendre coefficients";
  return os;
}

} // namespace som
