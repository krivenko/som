/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <string>

#include <triqs/gfs.hpp>

namespace som {

// All meshes used with input data containers
#define ALL_INPUT_MESHES (imtime)(imfreq)(legendre)

// Mesh traits
template <typename MeshType> struct mesh_traits;

template <> struct mesh_traits<triqs::gfs::imtime> {
  static constexpr int index = 0;
  inline static std::string name() { return "imaginary time"; }
};
template <> struct mesh_traits<triqs::gfs::imfreq> {
  static constexpr int index = 1;
  inline static std::string name() { return "imaginary frequency"; }
};
template <> struct mesh_traits<triqs::gfs::legendre> {
  static constexpr int index = 2;
  inline static std::string name() { return "Legendre"; }
};
template <> struct mesh_traits<triqs::gfs::refreq> {
  static constexpr int index = -1;
  inline static std::string name() { return "real frequency"; }
};

} // namespace som
