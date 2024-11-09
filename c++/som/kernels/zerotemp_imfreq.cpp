/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko
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

#include <numbers>

#include "zerotemp_imfreq.hpp"

namespace som {

using namespace triqs::gfs;

kernel<ZeroTemp, imfreq>::kernel(
    kernel<ZeroTemp, imfreq>::mesh_type const& mesh)
   : kernel_base(mesh.get_positive_freq().size())
   , mesh(mesh.get_positive_freq()) {}

// Apply to a rectangle
void kernel<ZeroTemp, imfreq>::apply(rectangle const& rect,
                                     result_view_type res) const {

  double e1 = rect.left();
  double e2 = rect.right();

  for(auto iw : mesh)
    res(iw.data_index()) =
        rect.height * std::log((dcomplex(iw) - e1) / (dcomplex(iw) - e2));
}

std::ostream& operator<<(std::ostream& os,
                         kernel<ZeroTemp, imfreq> const& kern) {
  os << R"(A(ϵ) -> G_{T=0}(iω), )";
  os << "Statistics = "
     << (kern.mesh.statistic() == Fermion ? "Fermion" : "Boson") << ", "
     << R"(Δω = )" << (2 * std::numbers::pi) / kern.mesh.beta() << ", "
     << kern.mesh.size() << " Matsubara frequencies";
  return os;
}

} // namespace som
