/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016 by I. Krivenko
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
#include <triqs/utility/numeric_ops.hpp>

#include "base.hpp"
#include "../spline.hpp"

#warning kernel<BosonCorr,imtime> is not implemented

namespace som {

using namespace triqs::arrays;
using namespace triqs::gfs;

// Kernel: bosonic correlator, imaginary time
template<> class kernel<BosonCorr,imtime> :
           public kernel_base<kernel<BosonCorr,imtime>, array<double,1>> {

public:

 using result_type = array<double,1>;
 using mesh_type = gf_mesh<imtime>;

 const double beta;          // Inverse temperature
 const mesh_type mesh;       // Matsubara time mesh

 kernel(mesh_type const& mesh) :
  kernel_base(mesh.size()), mesh(mesh), beta(mesh.domain().beta) {
  // TODO
 }

 // Apply to a rectangle
 void apply(rectangle const& rect, result_type & res) const {
  // TODO
 }

 friend std::ostream & operator<<(std::ostream & os, kernel const& kern) {
  os << "A(\\epsilon) -> \\chi(\\tau), ";
  os << "\\beta = " << kern.beta << ", " << kern.mesh.size() << " \\tau-points.";
  return os;
 }

};

}
