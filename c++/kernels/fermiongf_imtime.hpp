/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 by I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <cmath>
#include <triqs/gfs.hpp>

#include "base.hpp"
#include "../spline.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

// Kernel: fermionic GF, imaginary time
template<> class kernel<FermionGf,imtime> :
           public kernel_base<kernel<FermionGf,imtime>, array<double,1>> {

 // Tolerance levels for function evaluation
 static constexpr double tolerance = 1e-11;
 // Number of energy points for Lambda_tau_not0 interpolation
 static constexpr int n_spline_knots = 10001;

 // Integrated kernel \Lambda(\tau=0,\Omega)
 inline double Lambda_tau_0(double Omega) const {
  double x = beta*Omega;
  return (Omega < 0) ? (-std::log1p(std::exp(x))     )/beta :
                       (-std::log1p(std::exp(-x)) - x)/beta;
 }

 // Spline interpolation of the integrated kernel \Lambda(\tau!=0,\Omega)
 std::vector<regular_spline> Lambda_tau_not0;

public:

 using result_type = array<double,1>;
 using mesh_type = gf_mesh<imtime>;

 const double beta;          // Inverse temperature
 const gf_mesh<imtime> mesh; // Matsubara time mesh

 kernel(gf_mesh<imtime> const& mesh) : mesh(mesh), beta(mesh.x_max()) {

  Lambda_tau_not0.reserve(mesh.size()-2);
  for(int itau = 1; itau < mesh.size()-1; ++itau) {
   double alpha = mesh[itau] / beta;
   // Estimated limits of interpolation segment
   double Omega_min = (std::log(tolerance) + std::log(1-alpha)) / (1-alpha) / beta;
   double Omega_max = (std::log(tolerance) + std::log(alpha)) / (-alpha) / beta;
   double dOmega = (Omega_max - Omega_min) / (n_spline_knots - 1);

   // Integrated kernel \Lambda(\tau,\Omega=+\infty)
   double Lambda_inf = -M_PI/(beta * std::sin(M_PI * alpha));

   // Evaluate \Lambda(\tau,\Omega) to construct spline
   vector<double> values(n_spline_knots);
   for(int i = 0; i < n_spline_knots; ++i) {
    double Omega = Omega_min + i*dOmega;
    double val = 0;
    if(Omega > 0) {
     for(int n = 0;;++n) {
      double z = beta*(n + alpha);
      double t = (n % 2 ? 1 : -1) * std::exp(-Omega*z)/z;
      if(std::abs(t) < tolerance) break;
      val += t;
     }
     values[i] = Lambda_inf - val;
    } else if(Omega < 0) {
     for(int n = 0; ; ++n) {
      double z = beta*((n+1) - alpha);
      double t = (n % 2 ? 1 : -1) * std::exp(Omega*z)/z;
      if(std::abs(t) < tolerance) break;
      val += t;
     }
     values[i] = val;
    } else { // Omega == 0
     for(int n = 0; ; ++n) {
      double z = beta*((2*n+1) - alpha);
      double t = -beta / (z * (z + beta));
      if(std::abs(t) < tolerance) break;
      val += t;
     }
    values[i] = val;
    }
   }
   Lambda_tau_not0.push_back(regular_spline(Omega_min, Omega_max, values));
  }
 }

 void apply(rectangle const& rect, result_type & res) const {

  double e1 = rect.center - rect.width/2;
  double e2 = rect.center + rect.width/2;

  // (kernel * rect)(\tau = 0)
  res(0) = rect.height * (Lambda_tau_0(e2) - Lambda_tau_0(e1));
  // (kernel * rect)(0 < \tau < \beta)
  for(int itau = 1; itau < mesh.size()-1; ++itau) {
   auto const& l = Lambda_tau_not0[itau-1];
   res(itau) = rect.height * (l(e2) - l(e1));
  }
  // (kernel * rect)(\tau = \beta)
  res(mesh.size()-1) = rect.height * (Lambda_tau_0(-e1) - Lambda_tau_0(-e2));
 }

};

}
