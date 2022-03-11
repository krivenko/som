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

#include <iostream>
#include <string>

#include <som/kernels/mesh_traits.hpp>

#include "common.hxx"
#include "som_core.hpp"

namespace som {

using std::to_string;
using namespace triqs::gfs;

template <typename... GfOpts>
void check_input_gf(gf_const_view<GfOpts...> g, gf_const_view<GfOpts...> S) {
  if(g.mesh() != S.mesh() || g.target_shape() != S.target_shape())
    fatal_error(
        "input quantity and the error-bar function S must have equivalent "
        "structure");

  auto shape = g.target_shape();
  if(shape[0] != shape[1])
    fatal_error("matrix-valued input quantities must be square");
}

//////////////////////
// som_core::data_t //
//////////////////////

template <typename Mesh>
som_core::data_t::data_t(Mesh const& mesh, cache_index& ci)
   : rhs(input_data_t<Mesh>(mesh.size()))
   , error_bars(input_data_t<Mesh>(mesh.size()))
   , final_solution(ci) {}

template som_core::data_t::data_t(triqs::mesh::imtime const&, cache_index&);
template som_core::data_t::data_t(triqs::mesh::imfreq const&, cache_index&);
template som_core::data_t::data_t(triqs::mesh::legendre const&, cache_index&);

template <typename... GfOpts>
void som_core::data_t::init_input(int i, gf_const_view<GfOpts...> g,
                                  gf_const_view<GfOpts...> S, double norm_) {
  using mesh_t = typename gf_const_view<GfOpts...>::mesh_t;
  std::get<input_data_t<mesh_t>>(rhs)() = g.data()(range(), i, i);
  std::get<input_data_t<mesh_t>>(error_bars)() = S.data()(range(), i, i);
  if(norm_ <= 0)
    fatal_error("solution norm must be positive (got norm = " +
                std::to_string(norm_) + ")");
  norm = norm_;
}

////////////////////////////
// som_core: Constructors //
////////////////////////////

// Imaginary time
som_core::som_core(gf_const_view<imtime> g_tau, gf_const_view<imtime> S_tau,
                   observable_kind kind, vector<double> const& norms)
   : kind(kind)
   , mesh(g_tau.mesh())
   , data(g_tau.target_shape()[0], data_t(g_tau.mesh(), ci)) {

  if(is_stat_relevant(kind)) check_gf_stat(g_tau, observable_statistics(kind));

  check_input_gf(g_tau, S_tau);
  if(!is_gf_real(g_tau) || !is_gf_real(S_tau))
    fatal_error("imaginary time " + observable_name(kind) + " must be real");
  gf<imtime, matrix_real_valued> g_tau_real = real(g_tau),
                                 S_tau_real = real(S_tau);

  int gf_dim = g_tau.target_shape()[0];
  for(int i : range(gf_dim)) {
    data[i].init_input(i, make_const_view(g_tau_real),
                       make_const_view(S_tau_real),
                       norms.size() > 0 ? norms[i] : 1.0);
  }
}

// Imaginary frequency
som_core::som_core(gf_const_view<imfreq> g_iw, gf_const_view<imfreq> S_iw,
                   observable_kind kind, vector<double> const& norms)
   : kind(kind)
   , mesh(g_iw.mesh())
   , data(g_iw.target_shape()[0], data_t(positive_freq_view(g_iw).mesh(), ci)) {

  if(is_stat_relevant(kind)) check_gf_stat(g_iw, observable_statistics(kind));

  check_input_gf(g_iw, S_iw);
  if(!is_gf_real_in_tau(g_iw) || !is_gf_real_in_tau(S_iw))
    fatal_error("imaginary frequency " + observable_name(kind) +
                R"( must be real in \tau-domain)");
  auto g_iw_pos = positive_freq_view(g_iw);
  auto S_iw_pos = positive_freq_view(S_iw);
  check_input_gf(g_iw_pos, S_iw_pos);

  int gf_dim = g_iw_pos.target_shape()[0];
  for(int i : range(gf_dim)) {
    data[i].init_input(i, g_iw_pos, S_iw_pos,
                       norms.size() > 0 ? norms[i] : 1.0);
  }
}

// Legendre coefficients
som_core::som_core(gf_const_view<legendre> g_l, gf_const_view<legendre> S_l,
                   observable_kind kind, vector<double> const& norms)
   : kind(kind)
   , mesh(g_l.mesh())
   , data(g_l.target_shape()[0], data_t(g_l.mesh(), ci)) {

  if(is_stat_relevant(kind)) check_gf_stat(g_l, observable_statistics(kind));

  check_input_gf(g_l, S_l);
  if(!is_gf_real(g_l) || !is_gf_real(S_l))
    fatal_error("Legendre " + observable_name(kind) + " must be real");
  gf<legendre, matrix_real_valued> g_l_real = real(g_l), S_l_real = real(S_l);

  int gf_dim = g_l_real.target_shape()[0];
  for(int i : range(gf_dim)) {
    data[i].init_input(i, make_const_view(g_l_real), make_const_view(S_l_real),
                       norms.size() > 0 ? norms[i] : 1.0);
  }
}

} // namespace som
