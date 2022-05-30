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
#pragma once

#include <cmath>
#include <triqs/gfs.hpp>

#include "base.hpp"

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

// Kernel: bosonic correlator, Matsubara frequencies
template<> class kernel<BosonAutoCorr,imfreq> :
           public kernel_base<kernel<BosonAutoCorr,imfreq>, array<dcomplex,1>> {

public:

 using result_type = array<dcomplex,1>;
 using mesh_type = gf_mesh<imfreq>;
 constexpr static observable_kind kind = BosonAutoCorr;

 const double beta;          // Inverse temperature
 const mesh_type mesh; // Matsubara frequency mesh

 kernel(mesh_type const& mesh) :
  kernel_base(mesh.get_positive_freq().size()), mesh(mesh.get_positive_freq()),
  beta(mesh.domain().beta) {}

 // Apply to a rectangle
 void apply(rectangle const& rect, result_type & res) const {

  double e1 = rect.center - rect.width/2;
  double e2 = rect.center + rect.width/2;

  auto it = std::begin(mesh);
  res(it->index()) = 2 * rect.height * rect.width / M_PI; // \Omega = 0
  for(++it; it != std::end(mesh); ++it) {
   auto w = dcomplex(*it).imag();
   res(it->index()) = (2 * rect.height/M_PI) * (rect.width + w * (std::atan(e1/w) - std::atan(e2/w)));
  }
 }

 friend std::ostream & operator<<(std::ostream & os, kernel const& kern) {
  os << "A(\\epsilon) -> \\chi_{sym}(i\\Omega), ";
  os << "\\beta = " << kern.beta << ", " << kern.mesh.size() << " Matsubara frequencies.";
  return os;
 }

};

}
