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

#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/enum.hpp>

#include "../rectangle.hpp"
#include "../configuration.hpp"
#include "../config_update.hpp"

namespace som {

// Kinds of observables
#define ALL_OBSERVABLES (FermionGf)(BosonCorr)(BosonAutoCorr)(ZeroTemp)

enum observable_kind : unsigned int {BOOST_PP_SEQ_ENUM(ALL_OBSERVABLES)};
constexpr const unsigned int n_observable_kinds = BOOST_PP_SEQ_SIZE(ALL_OBSERVABLES);

// Is statistics defined for this observable?
inline bool is_stat_relevant(observable_kind kind) {
 return kind != ZeroTemp;
}

// Statistics of observables
inline triqs::gfs::statistic_enum observable_statistics(observable_kind kind) {
 switch(kind) {
  case FermionGf: return triqs::gfs::Fermion;
  case BosonCorr: return triqs::gfs::Boson;
  case BosonAutoCorr: return triqs::gfs::Boson;
  default: TRIQS_RUNTIME_ERROR <<
      "observable_statistics: Statistics is undefined for this observable kind";
 }
}

// Names of observables
inline std::string observable_name(observable_kind kind) {
 switch(kind) {
  case FermionGf: return "Green's functions";
  case BosonCorr: return "bosonic correlator";
  case BosonAutoCorr: return "bosonic autocorrelator";
  case ZeroTemp: return "correlator at zero temperature";
  default: TRIQS_RUNTIME_ERROR << "observable_name: Invalid observable kind";
 }
}

// Widest energy windows for observables
inline std::pair<double,double> max_energy_window(observable_kind kind) {
 switch(kind) {
  case FermionGf: return std::make_pair(-HUGE_VAL,HUGE_VAL);
  case BosonCorr: return std::make_pair(-HUGE_VAL,HUGE_VAL);
  case BosonAutoCorr: return std::make_pair(0,HUGE_VAL);
  case ZeroTemp: return std::make_pair(0,HUGE_VAL);
  default: TRIQS_RUNTIME_ERROR << "max_energy_window: Invalid observable kind";
 }
}

// Construct a real-frequency GF from a configuration
inline void back_transform(observable_kind kind,
                           configuration const& conf,
                           cache_index & ci,
                           triqs::gfs::gf_view<triqs::gfs::refreq,triqs::gfs::scalar_valued> g_w) {
 bool bosoncorr = kind == BosonCorr || kind == BosonAutoCorr;

 g_w() = 0;

 for(auto const& rect : conf) {
  for(auto e : g_w.mesh()) g_w.data()(e.index()) += rect.hilbert_transform(double(e), bosoncorr);

  if(kind == BosonAutoCorr) {
   // Add a reflected rectangle
   rectangle reflected_rect(-rect.center, rect.width, rect.height, ci);
   for(auto e : g_w.mesh()) g_w.data()(e.index()) += reflected_rect.hilbert_transform(double(e), true);
  }
 }

 if(bosoncorr) {
  g_w.data() *= -1.0/M_PI;
 }
}

// Compute the GF tail from a configuration
inline auto compute_tail(observable_kind kind,
                         configuration const& conf,
                         cache_index & ci,
                         int max_order) {
 assert(max_order >= 0);

 bool bosoncorr = kind == BosonCorr || kind == BosonAutoCorr;

 triqs::arrays::array<std::complex<double>, 1> tail(max_order + 1);
 tail() = .0;

 for(auto const& rect : conf) {
  tail += rect.tail_coefficients(0, max_order, bosoncorr);

  if(kind == BosonAutoCorr) {
   // Add a reflected rectangle
   rectangle reflected_rect(-rect.center, rect.width, rect.height, ci);
   tail += reflected_rect.tail_coefficients(0, max_order, true);
  }
 }

 if(bosoncorr) tail *= -1.0/M_PI;

 return tail;
}

// All meshes used with input data containers
#define ALL_INPUT_MESHES (imtime)(imfreq)(legendre)

// Mesh traits
template<typename MeshType> struct mesh_traits;

template<> struct mesh_traits<triqs::gfs::imtime> {
 static constexpr int index = 0;
 static std::string name() { return "imaginary time"; }
};
template<> struct mesh_traits<triqs::gfs::imfreq> {
 static constexpr int index = 1;
 static std::string name() { return "imaginary frequency"; }
};
template<> struct mesh_traits<triqs::gfs::legendre> {
 static constexpr int index = 2;
 static std::string name() { return "Legendre"; }
};
template<> struct mesh_traits<triqs::gfs::refreq> {
 static constexpr int index = -1;
 static std::string name() { return "real frequency"; }
};

// Integral kernels
template<observable_kind, typename Mesh> class kernel;

// Base class for all kernels (CRTP)
template<typename Derived, typename ResultType>
class kernel_base {

 // Precomputed LHS values
 mutable std::vector<ResultType> lhs_cache;

 inline void check_cache_size(cache_index const& ci) const {
   if(lhs_cache.size() < ci.size())
    lhs_cache.resize(ci.size(), ResultType(first_dim(lhs_cache[0])));
 }

public:

#define DERIVED static_cast<Derived const*>(this)

 kernel_base(int res_size) {
  lhs_cache.resize(CACHE_SIZE, ResultType(res_size));
 }

 inline ResultType const& operator()(rectangle const& rect) const {
  auto & ce = rect.ci[rect.cache_id];
  check_cache_size(rect.ci);
  auto & res = lhs_cache[rect.cache_id];
  // Do we have a precomputed LHS for rectangle r?
  if(!ce.valid) {
   DERIVED->apply(rect, res);
   ce.valid = true;
  }
  return res;
 }

 inline ResultType const& operator()(configuration const& c) const {
  auto & ce = c.ci[c.cache_id];
  check_cache_size(c.ci);
  auto & res = lhs_cache[c.cache_id];
  // Do we have a precomputed LHS for configuration c?
  if(!ce.valid) {
   res() = 0;
   for(auto const& r : c) res += operator()(r);
   ce.valid = true;
  }
  return res;
 }

 inline ResultType const& operator()(config_update const& cu) const {
  auto const& conf = cu.get_config();

  auto & ce = cu.ci[cu.cache_id];
  check_cache_size(cu.ci);
  auto & res = lhs_cache[cu.cache_id];
  // Do we have a precomputed LHS for config_update cu?
  if(!ce.valid) {
   // Cache entry holds the LHS for the *updated configuration*
   res = operator()(conf);

   auto rect_it = std::begin(cu.new_rects);
   for(int index : cu.changed_indices) {
    if(index == INT_MAX) { // add rectangle
     res += operator()(*rect_it);
     ++rect_it;
    } else if(index < 0) { // remove rectangle
     res -= operator()(conf[-index-1]);
    } else { // change rectangle
     res += operator()(*rect_it) - operator()(conf[index]);
     ++rect_it;
    }
   }
   ce.valid = true;
  }
  return res;
 }

#undef DERIVED

 // Copy a precomputed LHS value bound to one object to a cache entry bound to another object.
 // The destination cache entry will be marked as valid
 template<typename T1, typename T2>
 void cache_copy(T1 const& from, T2 const& to) const {
  check_cache_size(from.ci);
  lhs_cache[to.cache_id] = lhs_cache[from.cache_id];
  to.ci[to.cache_id].valid = true;
 }

 // Swap two LHS values each bound to a given object
 // The destination cache entry (second argument) will be marked as valid
 template<typename T1, typename T2>
 void cache_swap(T1 const& from, T2 const& to) const {
  check_cache_size(from.ci);
  std::swap(lhs_cache[to.cache_id], lhs_cache[from.cache_id]);
  to.ci[to.cache_id].valid = true;
 }

};

}
