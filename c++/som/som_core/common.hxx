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
#pragma once

#include <iostream>
#include <string>

#include <mpi/mpi.hpp>

#include <triqs/gfs.hpp>

#include <som/kernels/mesh_traits.hpp>
#include <som/kernels/observables.hpp>
#include <som/worker_parameters.hpp>

namespace som {

using std::to_string;
using namespace triqs::mesh;
using namespace triqs::gfs;

[[noreturn]] inline static void fatal_error(std::string const& message) {
  TRIQS_RUNTIME_ERROR << "som_core: " << message;
}

inline std::ostream& mpi_cout(mpi::communicator const& comm) {
  return (std::cout << "[Rank " << comm.rank() << "] ");
}

template <typename MeshType>
void check_gf_dim(gf_const_view<MeshType> g, std::size_t expected_dim) {
  auto shape = g.target_shape();
  if(shape[0] != expected_dim || shape[1] != expected_dim)
    fatal_error("expected a " + mesh_traits<MeshType>::name() +
                " Green's function with matrix dimensions " +
                to_string(expected_dim) + "x" + to_string(expected_dim));
}
template <typename MeshType>
void check_gf_dim(gf_view<MeshType> g, std::size_t expected_dim) {
  check_gf_dim(make_const_view(g), expected_dim);
}

template <typename MeshType>
void check_gf_stat(gf_const_view<MeshType> g, statistic_enum expected_stat) {
  if(g.mesh().statistic() != expected_stat)
    fatal_error("expected a " + mesh_traits<MeshType>::name() +
                " Green's function with " +
                (expected_stat == Fermion ? "fermionic" : "bosonic") +
                " statistics");
}
template <typename MeshType>
void check_gf_stat(gf_view<MeshType> g, statistic_enum expected_stat) {
  return check_gf_stat(make_const_view(g), expected_stat);
}

template <typename MeshType>
inline constexpr std::size_t kernel_id(observable_kind kind) {
  return std::size_t(kind) + n_observable_kinds * mesh_traits<MeshType>::index;
}

inline std::size_t kernel_id(observable_kind kind, mesh_variant_t const& mesh) {
  return std::size_t(kind) + n_observable_kinds * mesh.index();
}

#define FOR_EACH_KERNEL(F)                                                     \
  BOOST_PP_SEQ_FOR_EACH_PRODUCT(F, (ALL_OBSERVABLES)(ALL_INPUT_MESHES))

#define SELECT_KERNEL(F, F_NAME)                                               \
  switch(kernel_id(kind, mesh)) {                                              \
    FOR_EACH_KERNEL(F)                                                         \
    default:                                                                   \
      fatal_error("Unknown observable kind " + std::to_string(kind) +          \
                  " in " #F_NAME);                                             \
  }

} // namespace som
