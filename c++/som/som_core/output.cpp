/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2026 Igor Krivenko
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

#include <string>
#include <type_traits>
#include <vector>

#include <boost/preprocessor/seq/for_each.hpp>

#include <som/kernels/all.hpp>

#include "common.hxx"
#include "som_core.hpp"

namespace som {

using namespace triqs::gfs;

///////////////////
// fill_refreq() //
///////////////////

void fill_refreq(gf_view<refreq> g_w, som_core const& cont, bool with_binning) {
  auto kind = cont.get_observable_kind();
  auto gf_dim = cont.get_dim();
  check_gf_dim(g_w, gf_dim);
  for(auto i : range(gf_dim)) {
    back_transform(kind,
                   cont.get_solution(i),
                   slice_target_to_scalar(g_w, i, i),
                   with_binning,
                   cont.get_comm());
  }
}

void fill_refreq(gf_view<refreq> g_w,
                 observable_kind kind,
                 std::vector<configuration> const& solutions,
                 bool with_binning) {
  auto gf_dim = long(solutions.size());
  check_gf_dim(g_w, gf_dim);
  for(auto i : range(gf_dim)) {
    back_transform(
        kind, solutions[i], slice_target_to_scalar(g_w, i, i), with_binning);
  }
}

////////////////////
// compute_tail() //
////////////////////

array<dcomplex, 3> compute_tail(int max_order, som_core const& cont) {
  auto kind = cont.get_observable_kind();
  auto gf_dim = cont.get_dim();
  array<dcomplex, 3> tail(max_order + 1, gf_dim, gf_dim);
  for(auto i : range(gf_dim)) {
    tail(range::all, i, i) = som::compute_tail(
        kind, cont.get_solution(i), max_order, cont.get_comm());
  }
  return tail;
}

array<dcomplex, 3> compute_tail(int max_order,
                                observable_kind kind,
                                std::vector<configuration> const& solutions) {
  auto gf_dim = long(solutions.size());
  array<dcomplex, 3> tail(max_order + 1, gf_dim, gf_dim);
  for(auto i : range(gf_dim)) {
    tail(range::all, i, i) = som::compute_tail(kind, solutions[i], max_order);
  }
  return tail;
}

///////////////////
// reconstruct() //
///////////////////

void fill_data(gf_view<imtime> g_tau, long i, vector<double> const& data) {
  g_tau.data()(range::all, i, i) = data;
}

void fill_data(gf_view<imfreq> g_iw, long i, vector<dcomplex> const& data) {
  auto g_positive_freq = positive_freq_view(g_iw);
  g_positive_freq.data()(range::all, i, i) = data;
  g_iw = make_gf_from_real_gf(make_const_view(g_positive_freq));
}

void fill_data(gf_view<legendre> g_l, long i, vector<double> const& data) {
  g_l.data()(range::all, i, i) = data;
}

template <typename MeshType, typename Solutions>
void reconstruct_impl(gf_view<MeshType> g,
                      long gf_dim,
                      observable_kind kind,
                      Solutions const& sols) {
  check_gf_dim(g, gf_dim);
  if(is_stat_relevant(kind)) check_gf_stat(g, observable_statistics(kind));

  g() = 0;
#define FILL_DATA_CASE(r, d, ok)                                               \
  case ok: {                                                                   \
    kernel<ok, MeshType> kern(g.mesh());                                       \
    for(auto i : range(gf_dim)) {                                              \
      if constexpr(std::is_same_v<Solutions, som_core>)                        \
        fill_data(g, i, kern.apply_wo_caching(sols.get_solution(i)));          \
      else                                                                     \
        fill_data(g, i, kern.apply_wo_caching(sols[i]));                       \
    }                                                                          \
    return;                                                                    \
  }
  switch(kind) {
    BOOST_PP_SEQ_FOR_EACH(FILL_DATA_CASE, _, ALL_OBSERVABLES)
    default: fatal_error("unknown observable kind " + std::to_string(kind));
  }
#undef FILL_DATA_CASE
}

void reconstruct(gf_view<imtime> g, som_core const& cont) {
  reconstruct_impl(g, cont.get_dim(), cont.get_observable_kind(), cont);
}
void reconstruct(gf_view<imfreq> g, som_core const& cont) {
  reconstruct_impl(g, cont.get_dim(), cont.get_observable_kind(), cont);
}
void reconstruct(gf_view<legendre> g, som_core const& cont) {
  reconstruct_impl(g, cont.get_dim(), cont.get_observable_kind(), cont);
}

void reconstruct(gf_view<imtime> g,
                 observable_kind kind,
                 std::vector<configuration> const& solutions) {
  reconstruct_impl(g, long(solutions.size()), kind, solutions);
}
void reconstruct(gf_view<imfreq> g,
                 observable_kind kind,
                 std::vector<configuration> const& solutions) {
  reconstruct_impl(g, long(solutions.size()), kind, solutions);
}
void reconstruct(gf_view<legendre> g,
                 observable_kind kind,
                 std::vector<configuration> const& solutions) {
  reconstruct_impl(g, long(solutions.size()), kind, solutions);
}

} // namespace som
