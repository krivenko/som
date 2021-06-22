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

#include <iostream>

#include <triqs/gfs.hpp>

#include <som/kernels/mesh_traits.hpp>

namespace som {

using std::to_string;
using namespace triqs::gfs;

inline void fatal_error(std::string const& message) {
  TRIQS_RUNTIME_ERROR << "som_core: " << message;
}

inline void warning(std::string const& message) {
  std::cout << "WARNING: " << message << std::endl;
}

template <typename MeshType>
void check_gf_dim(gf_const_view<MeshType> g, int expected_dim) {
  auto shape = g.target_shape();
  if(shape[0] != expected_dim || shape[1] != expected_dim)
    fatal_error("expected a " + mesh_traits<MeshType>::name() +
                " Green's function with matrix dimensions " +
                to_string(expected_dim) + "x" + to_string(expected_dim));
}
template <typename MeshType>
void check_gf_dim(gf_view<MeshType> g, int expected_dim) {
  check_gf_dim(make_const_view(g), expected_dim);
}

template <typename MeshType>
void check_gf_stat(gf_const_view<MeshType> g, statistic_enum expected_stat) {
  if(g.domain().statistic != expected_stat)
    fatal_error("expected a " + mesh_traits<MeshType>::name() +
                " Green's function with " +
                (expected_stat == Fermion ? "fermionic" : "bosonic") +
                " statistics");
}
template <typename MeshType>
void check_gf_stat(gf_view<MeshType> g, statistic_enum expected_stat) {
  return check_gf_stat(make_const_view(g), expected_stat);
}

} // namespace som
