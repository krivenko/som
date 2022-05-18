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

#include "fermiongfsymm_imfreq.hpp"

namespace som {

using namespace triqs::gfs;

kernel<FermionGfSymm, imfreq>::kernel(
    kernel<FermionGfSymm, imfreq>::mesh_type const& mesh)
   : kernel_base(mesh.get_positive_freq().size())
   , beta(mesh.domain().beta)
   , mesh(mesh.get_positive_freq()) {}

// Apply to a rectangle
void kernel<FermionGfSymm, imfreq>::apply(rectangle const& rect,
                                          result_view_type res) const {

  double e1 = rect.left();
  double e2 = rect.right();

  for(auto iw : mesh) {
    auto w = dcomplex(iw).imag();
    res(iw.linear_index()) =
        1i * rect.height * (std::atan(e1 / w) - std::atan(e2 / w));
  }
}

std::ostream& operator<<(std::ostream& os,
                         kernel<FermionGfSymm, imfreq> const& kern) {
  os << R"(A(ϵ) -> G_{sym}(iω), )";
  os << R"(β = )" << kern.beta << ", " << kern.mesh.size()
     << " Matsubara frequencies";
  return os;
}

} // namespace som
