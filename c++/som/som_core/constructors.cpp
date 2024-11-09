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

#include <cmath>
#include <iostream>
#include <string>
#include <type_traits>

#include <triqs/utility/is_complex.hpp>

#include <som/kernels/mesh_traits.hpp>

#include "common.hxx"
#include "som_core.hpp"

namespace som {

using std::to_string;
using namespace triqs::mesh;

template <typename A> bool isfinite(A const& arr) {
  using value_type = std::remove_const_t<typename A::value_type>;
  for(auto x : arr) {
    if constexpr(triqs::is_complex<value_type>::value) {
      if((!std::isfinite(x.real())) || (!std::isfinite(x.imag()))) return false;
    } else {
      if(!std::isfinite(x)) return false;
    }
  }
  return true;
}

template <typename Mesh, typename... GfOpts>
void check_input_gf_and_error_bars(gf_const_view<Mesh, GfOpts...> g,
                                   gf_const_view<Mesh, GfOpts...> error_bars) {
  if(g.mesh() != error_bars.mesh() ||
     g.target_shape() != error_bars.target_shape())
    fatal_error(
        "input quantity and the error bars function must have equivalent "
        "structures");

  auto shape = g.target_shape();
  if(shape[0] != shape[1])
    fatal_error("matrix-valued input quantities must be square");

  if(!isfinite(g.data()))
    fatal_error("input quantity contains infinite or NaN values");
  if(!isfinite(error_bars.data()))
    fatal_error("error bars contain infinite or NaN values");
}

template <typename Mesh, typename... GfOpts>
void check_input_gf_and_cov_matrices(
    gf_const_view<Mesh, GfOpts...> g,
    som_core::cov_matrices_gf_view_type<Mesh> cov_matrices) {
  auto shape = g.target_shape();
  if(shape[0] != shape[1])
    fatal_error("matrix-valued input quantities must be square");

  if(cov_matrices.target_shape()[0] != shape[0])
    fatal_error("Exactly " + std::to_string(shape[0]) +
                " covariance matrices must be provided");

  if(std::get<0>(cov_matrices.mesh()) != g.mesh() ||
     std::get<1>(cov_matrices.mesh()) != g.mesh())
    fatal_error(
        "Covariance matrices are defined on a mesh different from that of the "
        "RHS");

  if(!isfinite(g.data()))
    fatal_error("input quantity contains infinite or NaN values");
  if(!isfinite(cov_matrices.data()))
    fatal_error("covariance matrices contain infinite or NaN values");
}

void check_norms(nda::vector<double> const& norms, long gf_dim) {
  if(norms.size() != gf_dim)
    fatal_error("The 'norms' list must have exactly " + std::to_string(gf_dim) +
                " elements (got " + std::to_string(norms.size()) +
                " elements)");
}

void check_filtering_levels(nda::vector<double> const& filtering_levels,
                            long gf_dim) {
  if(filtering_levels.size() > 0 && filtering_levels.size() != gf_dim)
    fatal_error(
        "The 'filtering_levels' list must either be empty or have exactly " +
        std::to_string(gf_dim) + " elements (got " +
        std::to_string(filtering_levels.size()) + " elements)");
}

//////////////////////
// som_core::data_t //
//////////////////////

template <typename Mesh, typename TargetOpt>
som_core::data_t::data_t(Mesh const& mesh,
                         cache_index& ci,
                         triqs::gfs::gf_const_view<Mesh, TargetOpt> g,
                         triqs::gfs::gf_const_view<Mesh, TargetOpt> error_bars,
                         double norm)
   : rhs(input_data_t<Mesh>(mesh.size()))
   , errors(input_data_t<Mesh>(mesh.size()))
   , norm(norm)
   , final_solution(ci) {

  for(auto sigma : error_bars.data()) {
    if(std::real(sigma) <= 0) fatal_error("All error bars must be positive");
  }

  std::get<input_data_t<Mesh>>(rhs)() = g.data();
  std::get<input_data_t<Mesh>>(errors)() = error_bars.data();
  if(norm <= 0)
    fatal_error("solution norm must be positive (got norm = " +
                std::to_string(norm) + ")");
}

template <typename Mesh, typename TargetOpt, typename CovMatrixTargetOpt>
som_core::data_t::data_t(
    Mesh const& mesh,
    cache_index& ci,
    triqs::gfs::gf_const_view<Mesh, TargetOpt> g,
    triqs::gfs::gf_const_view<prod<Mesh, Mesh>, CovMatrixTargetOpt> cov_matrix,
    double norm,
    double filtering_level)
   : rhs(input_data_t<Mesh>(mesh.size()))
   , errors(cov_matrix_t<Mesh>((long)mesh.size(), (long)mesh.size()))
   , filtering_level(filtering_level)
   , norm(norm)
   , final_solution(ci) {

  if(dagger(cov_matrix.data()) != cov_matrix.data())
    fatal_error("Covariance matrix is not Hermitian");

  std::get<input_data_t<Mesh>>(rhs)() = g.data();
  std::get<cov_matrix_t<Mesh>>(errors)() = cov_matrix.data();
  if(norm <= 0)
    fatal_error("solution norm must be positive (got norm = " +
                std::to_string(norm) + ")");
}

////////////////////////////
// som_core: Constructors //
////////////////////////////

//
// Imaginary time
//

som_core::som_core(triqs::gfs::gf_const_view<imtime> g_tau,
                   triqs::gfs::gf_const_view<imtime> error_bars_tau,
                   observable_kind kind,
                   nda::vector<double> const& norms)
   : kind(kind), mesh(g_tau.mesh()) {

  if(is_stat_relevant(kind)) check_gf_stat(g_tau, observable_statistics(kind));

  check_input_gf_and_error_bars(g_tau, error_bars_tau);
  if(!is_gf_real(g_tau))
    fatal_error("imaginary time " + observable_name(kind) + " must be real");
  if(!is_gf_real(error_bars_tau)) fatal_error("error bars must be real");
  auto const g_tau_real = real(g_tau);
  auto const error_bars_tau_real = real(error_bars_tau);

  auto gf_dim = g_tau.target_shape()[0];

  check_norms(norms, gf_dim);

  data.reserve(gf_dim);
  for(auto n : range(gf_dim)) {
    // cppcheck-suppress useStlAlgorithm
    data.emplace_back(g_tau.mesh(),
                      ci,
                      slice_target_to_scalar(g_tau_real, n, n),
                      slice_target_to_scalar(error_bars_tau_real, n, n),
                      norms[n]);
  }
}

som_core::som_core(triqs::gfs::gf_const_view<imtime> g_tau,
                   cov_matrices_gf_view_type<imtime> cov_matrices_tau,
                   observable_kind kind,
                   nda::vector<double> const& norms,
                   nda::vector<double> const& filtering_levels)
   : kind(kind), mesh(g_tau.mesh()) {

  if(is_stat_relevant(kind)) check_gf_stat(g_tau, observable_statistics(kind));

  check_input_gf_and_cov_matrices(g_tau, cov_matrices_tau);
  if(!is_gf_real(g_tau))
    fatal_error("imaginary time " + observable_name(kind) + " must be real");
  auto const g_tau_real = real(g_tau);

  auto gf_dim = g_tau.target_shape()[0];

  check_norms(norms, gf_dim);
  check_filtering_levels(filtering_levels, gf_dim);

  data.reserve(gf_dim);
  for(auto n : range(gf_dim)) {
    auto const cov_matrix_real =
        real(slice_target_to_scalar(cov_matrices_tau, n));
    // cppcheck-suppress useStlAlgorithm
    data.emplace_back(g_tau.mesh(),
                      ci,
                      slice_target_to_scalar(g_tau_real, n, n),
                      cov_matrix_real(),
                      norms[n],
                      filtering_levels.size() > 0 ? filtering_levels[n] : 0.0);
  }
}

//
// Imaginary frequency
//

som_core::som_core(gf_const_view<imfreq> g_iw,
                   gf_const_view<imfreq> error_bars_iw,
                   observable_kind kind,
                   vector<double> const& norms)
   : kind(kind), mesh(g_iw.mesh()) {

  if(is_stat_relevant(kind)) check_gf_stat(g_iw, observable_statistics(kind));

  check_input_gf_and_error_bars(g_iw, error_bars_iw);
  if(!is_gf_real_in_tau(g_iw))
    fatal_error("imaginary frequency " + observable_name(kind) +
                " must be real in τ-domain");
  if(!is_gf_real(error_bars_iw)) fatal_error("error bars must be real");
  auto const g_iw_pos = positive_freq_view(g_iw);
  auto const error_bars_iw_pos = positive_freq_view(error_bars_iw);

  auto gf_dim = g_iw_pos.target_shape()[0];

  check_norms(norms, gf_dim);

  data.reserve(gf_dim);
  for(auto n : range(gf_dim)) {
    // cppcheck-suppress useStlAlgorithm
    data.emplace_back(g_iw_pos.mesh(),
                      ci,
                      slice_target_to_scalar(g_iw_pos, n, n),
                      slice_target_to_scalar(error_bars_iw_pos, n, n),
                      norms[n]);
  }
}

som_core::som_core(triqs::gfs::gf_const_view<imfreq> g_iw,
                   cov_matrices_gf_view_type<imfreq> cov_matrices_iw,
                   observable_kind kind,
                   nda::vector<double> const& norms,
                   nda::vector<double> const& filtering_levels)
   : kind(kind), mesh(g_iw.mesh()) {

  if(is_stat_relevant(kind)) check_gf_stat(g_iw, observable_statistics(kind));

  check_input_gf_and_cov_matrices(g_iw, cov_matrices_iw);
  if(!is_gf_real_in_tau(g_iw))
    fatal_error("imaginary frequency " + observable_name(kind) +
                " must be real in τ-domain");
  auto const g_iw_pos = positive_freq_view(g_iw);

  auto gf_dim = g_iw.target_shape()[0];

  check_norms(norms, gf_dim);
  check_filtering_levels(filtering_levels, gf_dim);

  auto make_cov_matrix_pos = [&](auto const& cov_matrix_iw) {
    using namespace triqs::gfs;
    auto cov_matrix_pos = gf<prod<imfreq, imfreq>, scalar_valued>{
        {g_iw_pos.mesh(), g_iw_pos.mesh()}};
    triqs::clef::placeholder<0> iw1;
    triqs::clef::placeholder<1> iw2;
    cov_matrix_pos(iw1, iw2) << cov_matrix_iw(iw1, iw2);
    return cov_matrix_pos;
  };

  data.reserve(gf_dim);
  for(auto n : range(gf_dim)) {
    auto const cov_matrix_pos =
        make_cov_matrix_pos(slice_target_to_scalar(cov_matrices_iw, n));
    // cppcheck-suppress useStlAlgorithm
    data.emplace_back(g_iw_pos.mesh(),
                      ci,
                      slice_target_to_scalar(g_iw_pos, n, n),
                      cov_matrix_pos(),
                      norms[n],
                      filtering_levels.size() > 0 ? filtering_levels[n] : 0.0);
  }
}

//
// Legendre coefficients
//

som_core::som_core(triqs::gfs::gf_const_view<legendre> g_l,
                   triqs::gfs::gf_const_view<legendre> error_bars_l,
                   observable_kind kind,
                   nda::vector<double> const& norms)
   : kind(kind), mesh(g_l.mesh()) {

  if(is_stat_relevant(kind)) check_gf_stat(g_l, observable_statistics(kind));

  check_input_gf_and_error_bars(g_l, error_bars_l);
  if(!is_gf_real(g_l))
    fatal_error("Legendre " + observable_name(kind) + " must be real");
  if(!is_gf_real(error_bars_l)) fatal_error("error bars must be real");
  auto const g_l_real = real(g_l);
  auto const error_bars_l_real = real(error_bars_l);

  auto gf_dim = g_l.target_shape()[0];

  check_norms(norms, gf_dim);

  data.reserve(gf_dim);
  for(auto n : range(gf_dim)) {
    // cppcheck-suppress useStlAlgorithm
    data.emplace_back(g_l.mesh(),
                      ci,
                      slice_target_to_scalar(g_l_real, n, n),
                      slice_target_to_scalar(error_bars_l_real, n, n),
                      norms[n]);
  }
}

som_core::som_core(triqs::gfs::gf_const_view<legendre> g_l,
                   cov_matrices_gf_view_type<legendre> cov_matrices_l,
                   observable_kind kind,
                   nda::vector<double> const& norms,
                   nda::vector<double> const& filtering_levels)
   : kind(kind), mesh(g_l.mesh()) {

  if(is_stat_relevant(kind)) check_gf_stat(g_l, observable_statistics(kind));

  check_input_gf_and_cov_matrices(g_l, cov_matrices_l);
  if(!is_gf_real(g_l))
    fatal_error("Legendre " + observable_name(kind) + " must be real");
  auto const g_l_real = real(g_l);

  auto gf_dim = g_l.target_shape()[0];

  check_norms(norms, gf_dim);
  check_filtering_levels(filtering_levels, gf_dim);

  data.reserve(gf_dim);
  for(auto n : range(gf_dim)) {
    auto const cov_matrix_real =
        real(slice_target_to_scalar(cov_matrices_l, n));
    // cppcheck-suppress useStlAlgorithm
    data.emplace_back(g_l.mesh(),
                      ci,
                      slice_target_to_scalar(g_l_real, n, n),
                      cov_matrix_real(),
                      norms[n],
                      filtering_levels.size() > 0 ? filtering_levels[n] : 0.0);
  }
}

} // namespace som
