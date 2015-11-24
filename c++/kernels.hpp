#pragma once

#include <boost/tr1/cmath.hpp>

#include <triqs/arrays/vector.hpp>
#include <triqs/gfs/meshes/matsubara_time.hpp>

#include "configuration.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

// Kinds of observables
enum observable_kind {greens_function, susceptibility, conductivity};

// Integral kernels
template<typename Mesh, statistic_enum Statistics, observable_kind Kind = greens_function>
class kernel;

template<>
class kernel<imtime,Fermion,greens_function> {

 double beta;          // Inverse temperature
 gf_mesh<imtime> mesh; // Matsubara time mesh

 kernel(gf_mesh<imtime> const& mesh) : mesh(mesh), beta(mesh.x_max()) {}

 vector<double> operator()(rectangle const& rect) {
  double e1 = rect.center - rect.width/2;
  double e2 = rect.center + rect.width/2;
  vector<double> res(mesh.size());

  auto mesh_it = std::begin(mesh);
  // (kernel * rect)(tau = 0):
  // -h *( (e2 + (1/beta) * 
  ++mesh_it;

  for(auto const& tau : mesh) {
  }

  return res;
 }

};

}
