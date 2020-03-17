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

#include <complex>
#include <iostream>

#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>

#include "base.hpp"

namespace som {

// Kernel: bosonic correlator, Matsubara frequencies
template <>
class kernel<BosonCorr, triqs::gfs::imfreq>
   : public kernel_base<kernel<BosonCorr, triqs::gfs::imfreq>,
                        triqs::arrays::array<std::complex<double>, 1>> {

public:
  using result_type = triqs::arrays::array<std::complex<double>, 1>;
  using mesh_type = triqs::gfs::gf_mesh<triqs::gfs::imfreq>;
  constexpr static observable_kind kind = BosonCorr;

  const double beta;    // Inverse temperature
  const mesh_type mesh; // Matsubara frequency mesh

  explicit kernel(mesh_type const& mesh);

  // Apply to a rectangle
  void apply(rectangle const& rect, result_type& res) const;

  friend std::ostream& operator<<(std::ostream& os, kernel const& kern);
};

} // namespace som
