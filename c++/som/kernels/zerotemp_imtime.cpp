/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2025 Igor Krivenko
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

#include "zerotemp_imtime.hpp"

namespace som {

using namespace triqs::gfs;

//////////////////////////////
// kernel<ZeroTemp, imtime> //
//////////////////////////////

kernel<ZeroTemp, imtime>::kernel(mesh_type const& mesh)
   : kernel_base(mesh.size()), mesh(mesh) {}

void kernel<ZeroTemp, imtime>::apply(rectangle const& rect,
                                     result_view_type res) const {

  double e1 = rect.left();
  double e2 = rect.right();

  auto it = mesh.begin();
  res(0) = -rect.height * rect.width;
  for(++it; it != mesh.end(); ++it) {
    auto tau = double(*it);
    res((*it).data_index()) =
        rect.height * (std::exp(-tau * e2) - std::exp(-tau * e1)) / tau;
  }
}

std::ostream& operator<<(std::ostream& os,
                         kernel<ZeroTemp, imtime> const& kern) {
  os << R"(A(ϵ) -> G_{T=0}(τ), )";
  os << "Statistics = "
     << (kern.mesh.statistic() == Fermion ? "Fermion" : "Boson") << ", "
     << R"(τ_{max} = )" << kern.mesh.beta() << ", " << kern.mesh.size()
     << R"( τ-points)";
  return os;
}

} // namespace som
