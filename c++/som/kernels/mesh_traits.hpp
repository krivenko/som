/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include <variant>

#include <triqs/mesh.hpp>

namespace som {

// All meshes used with input data containers
#define ALL_INPUT_MESHES (imtime)(imfreq)(legendre)

// Mesh traits
template <typename MeshType> struct mesh_traits;

template <> struct mesh_traits<triqs::mesh::imtime> {
  static constexpr std::size_t index = 0;
  inline static std::string name() { return "imaginary time"; }
  static constexpr bool is_complex_input_data = false;
};
template <> struct mesh_traits<triqs::mesh::imfreq> {
  static constexpr std::size_t index = 1;
  inline static std::string name() { return "imaginary frequency"; }
  static constexpr bool is_complex_input_data = true;
};
template <> struct mesh_traits<triqs::mesh::legendre> {
  static constexpr std::size_t index = 2;
  inline static std::string name() { return "Legendre"; }
  static constexpr bool is_complex_input_data = false;
};
template <> struct mesh_traits<triqs::mesh::refreq> {
  static constexpr std::size_t index = -1;
  inline static std::string name() { return "real frequency"; }
};

using mesh_variant_t = std::
    variant<triqs::mesh::imtime, triqs::mesh::imfreq, triqs::mesh::legendre>;

} // namespace som
