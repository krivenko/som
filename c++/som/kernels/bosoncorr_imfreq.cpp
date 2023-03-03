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

#include "bosoncorr_imfreq.hpp"

namespace som {

using namespace triqs::gfs;

kernel<BosonCorr, imfreq>::kernel(
    kernel<BosonCorr, imfreq>::mesh_type const& mesh)
   : kernel_base(mesh.get_positive_freq().size())
   , beta(mesh.domain().beta)
   , mesh(mesh.get_positive_freq()) {}

// Apply to a rectangle
void kernel<BosonCorr, imfreq>::apply(rectangle const& rect,
                                      result_view_type res) const {

  double e1 = rect.left();
  double e2 = rect.right();

  auto it = std::begin(mesh);
  res((*it).linear_index()) = rect.height * rect.width / M_PI; // \Omega = 0
  for(++it; it != std::end(mesh); ++it) {
    auto iw = dcomplex(*it);
    res((*it).linear_index()) =
        (rect.height / M_PI) *
        (rect.width + iw * std::log((iw - e2) / (iw - e1)));
  }
}

std::ostream& operator<<(std::ostream& os,
                         kernel<BosonCorr, imfreq> const& kern) {
  os << R"(A(ϵ) -> χ(iΩ), )";
  os << R"(β = )" << kern.beta << ", " << kern.mesh.size()
     << " Matsubara frequencies";
  return os;
}

} // namespace som
