/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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

// Kernel: zero temperature GF, imaginary time
template<> class kernel<ZeroTemp,imtime> :
           public kernel_base<kernel<ZeroTemp,imtime>, array<double,1>> {

public:

 using result_type = array<double,1>;
 using mesh_type = gf_mesh<imtime>;
 constexpr static observable_kind kind = ZeroTemp;

 const double beta = HUGE_VAL; // Inverse temperature (infinity)
 const mesh_type mesh;         // Matsubara time mesh

 kernel(mesh_type const& mesh) :
  kernel_base(mesh.size()), mesh(mesh) {}

 // Apply to a rectangle
 void apply(rectangle const& rect, result_type & res) const {

  double e1 = rect.center - rect.width/2;
  double e2 = rect.center + rect.width/2;

  auto it = mesh.begin();
  res(0) = -rect.height * rect.width;
  for(++it; it != mesh.end(); ++it) {
   auto tau = double(*it);
   res(it->index()) = rect.height * (std::exp(-tau*e2) - std::exp(-tau*e1)) / tau;
  }
 }

 friend std::ostream & operator<<(std::ostream & os, kernel const& kern) {
  os << "A(\\epsilon) -> G_{T=0}(\\tau), ";
  os << "Statistics = " << (kern.mesh.domain().statistic == Fermion ? "Fermion" : "Boson") << ", "
     << "\\tau_{max} = " << kern.mesh.domain().beta << ", "
     << kern.mesh.size() << " Matsubara frequencies.";
  return os;
 }

};

}
