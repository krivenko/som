#pragma once

#include <cmath>
#include <triqs/arrays/vector.hpp>
#include <triqs/gfs.hpp>

#include "configuration.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

// Kinds of observables
enum observable_kind {FermionicGf, Susceptibility, Conductivity};

// Integral kernels
template<observable_kind Kind, typename Mesh> class kernel;

// Kernel: fermionic GF, imaginary time
template<> class kernel<FermionicGf,imtime> {

 double beta;          // Inverse temperature
 gf_mesh<imtime> mesh; // Matsubara time mesh

 // Tolerance levels for function evaluation
 static constexpr double tolerance = 1e-10;

 // Integrated kernel \Lambda(\tau!=0,0)
 vector<double> Lambda_0;
 // Integrated kernel \Lambda(\tau!=0,+\infty)
 vector<double> Lambda_inf;
 // Integrated kernel \Lambda(\tau=0,\Omega)
 double Lambda_tau_0(double Omega) const {
  return (Omega < 0) ? -std::log(1.0 + std::exp(beta*Omega))/beta :
                       -std::log(1.0 + std::exp(-beta*Omega))/beta -Omega;
 }
 // Integrated kernel \Lambda(\tau=\beta,\Omega)
 // In this case the definition involves an integration from \Omega to +\infty
 double Lambda_tau_beta(double Omega) const {
  return (Omega < 0) ? -std::log(1.0 + std::exp(beta*Omega)) +Omega :
                       -std::log(1.0 + std::exp(-beta*Omega));
 }
 // Integrated kernel \Lambda(\tau!=0,\Omega)
 double Lambda_tau_not0(int itau, double Omega) const {
  if(Omega > 0) {
   double s = 0;
   for(long n = 0;;++n) {
    double z = beta*n + mesh[itau];
    double t = (n % 2 ? 1 : -1) * std::exp(-Omega*z)/z;
    if(std::abs(t) < tolerance) break;
    s += t;
   }
   return Lambda_inf[itau] - s;
  } else if(Omega < 0) {
   double s = 0;
   for(long n = 0; ; ++n) {
    double z = beta*(n+1) - mesh[itau];
    double t = (n % 2 ? 1 : -1) * std::exp(Omega*z)/z;
    if(std::abs(t) < tolerance) break;
    s += t;
   }
   return s;
  } else
   return Lambda_0[itau];
 }

public:

 kernel(gf_mesh<imtime> const& mesh) :
  mesh(mesh), beta(mesh.x_max()),
  Lambda_0(mesh.size()), Lambda_inf(mesh.size()) {
  // Fill Lambda_0
  for(int itau = 1; itau < mesh.size(); ++itau) {
   double s = 0;
   for(long n = 0; ; ++n) {
    double x = beta*(2*n+1) - mesh[itau];
    double t = -beta / (x * (x + beta));
    if(std::abs(t) < tolerance) break;
    s += t;
   }
   Lambda_0[itau] = s;
  }
  // Fill Lambda_inf
  for(int itau = 1; itau < mesh.size(); ++itau) {
   Lambda_inf[itau] = -M_PI/(beta * std::sin(M_PI * mesh[itau]/beta));
  }
 }

 vector<double> operator()(rectangle const& rect) const {
  double e1 = rect.center - rect.width/2;
  double e2 = rect.center + rect.width/2;
  vector<double> res(mesh.size());

  // (kernel * rect)(\tau = 0)
  res[0] = rect.height * (Lambda_tau_0(e2) - Lambda_tau_0(e1));
  // (kernel * rect)(0 < \tau < \beta)
  for(int itau = 1; itau < mesh.size()-1; ++itau) {
   res[itau] = rect.height * (Lambda_tau_not0(itau,e2) - Lambda_tau_not0(itau,e1));
  }
  res[mesh.size()-1] = rect.height * (Lambda_tau_beta(e1) - Lambda_tau_beta(e2));

  return res;
 }

 vector<double> operator()(configuration const& c) {
  vector<double> res(mesh.size());
  for(auto const& r : c.rects) res += operator()(r);
  return res;
 }

};

}
