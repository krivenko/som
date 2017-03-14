/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2017 by I. Krivenko
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

#include <cmath>
#include <triqs/gfs.hpp>

#include "base.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

// Kernel: zero temperature GF, Matsubara frequencies
template<> class kernel<ZeroTemp,imfreq> :
           public kernel_base<kernel<ZeroTemp,imfreq>, array<dcomplex,1>> {

public:

 using result_type = array<dcomplex,1>;
 using mesh_type = gf_mesh<imfreq>;
 constexpr static observable_kind kind = ZeroTemp;

 const double beta = HUGE_VAL; // Inverse temperature (infinity)
 const mesh_type mesh;         // Matsubara frequency mesh

 kernel(mesh_type const& mesh) :
  kernel_base(mesh.get_positive_freq().size()), mesh(mesh.get_positive_freq()) {}

 // Apply to a rectangle
 void apply(rectangle const& rect, result_type & res) const {

  double e1 = rect.center - rect.width/2;
  double e2 = rect.center + rect.width/2;

  for(auto iw : mesh)
   res(iw.index()) = rect.height * std::log((dcomplex(iw) - e1) / (dcomplex(iw) - e2));
 }

 friend std::ostream & operator<<(std::ostream & os, kernel const& kern) {
  os << "A(\\epsilon) -> G_{T=0}(i\\omega), ";
  os << "Statistics = " << (kern.mesh.domain().statistic == Fermion ? "Fermion" : "Boson") << ", "
     << "\\Delta\\omega = " << (2*M_PI)/kern.mesh.domain().beta << ", "
     << kern.mesh.size() << " Matsubara frequencies.";
  return os;
 }

};

}
