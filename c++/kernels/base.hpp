/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 by I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <string>
#include <vector>
#include <iostream>

#include "../rectangle.hpp"
#include "../configuration.hpp"
#include "../config_update.hpp"

namespace som {

// Kinds of observables
enum observable_kind {FermionGf, Susceptibility, Conductivity};

// Statistics of observables
triqs::gfs::statistic_enum observable_statistics(observable_kind kind) {
 switch(kind) {
  case FermionGf: return triqs::gfs::Fermion;
  case Susceptibility: return triqs::gfs::Boson;
  case Conductivity: return triqs::gfs::Boson;
 }
}

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

public:

#define DERIVED static_cast<Derived const*>(this)

 kernel_base(int res_size) {
  lhs_cache.resize(CACHE_SIZE, ResultType(res_size));
 };

 inline ResultType const& operator()(rectangle const& rect) const {
  auto & ce = *rect.cache_entry;
  auto & res = lhs_cache[ce.id];
  // Do we have a precomputed LHS for rectangle r?
  if(!ce.valid) {
   DERIVED->apply(rect, res);
   ce.valid = true;
  }
  return res;
 }

 inline ResultType const& operator()(configuration const& c) const {
  auto ce = c.cache_entry;
  auto & res = lhs_cache[ce->id];
  // Do we have a precomputed LHS for configuration c?
  if(!ce->valid) {
   res() = 0;
   for(auto const& r : c) res += operator()(r);
   ce->valid = true;
  }
  return res;
 }

 inline ResultType const& operator()(config_update const& cu) const {
  auto const& conf = cu.get_config();

  auto ce = cu.cache_entry;
  auto & res = lhs_cache[ce->id];
  // Do we have a precomputed LHS for config_update cu?
  if(!ce->valid) {
   // Cache entry holds the LHS for the *updated configuration*
   res = operator()(conf);

   auto rect_it = std::cbegin(cu.new_rects);
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
   ce->valid = true;
  }
  return res;
 }

#undef DERIVED

 // Copy a precomputed LHS value bound to one object to a cache entry bound to another object.
 // The destination cache entry will be marked as valid
 template<typename T1, typename T2>
 void cache_copy(T1 const& from, T2 const& to) const {
  lhs_cache[to.cache_entry->id] = lhs_cache[from.cache_entry->id];
  to.cache_entry->valid = true;
 }

 // Swap two LHS values each bound to a given object
 // The destination cache entry (second argument) will be marked as valid
 template<typename T1, typename T2>
 void cache_swap(T1 const& from, T2 const& to) const {
  std::swap(lhs_cache[to.cache_entry->id], lhs_cache[from.cache_entry->id]);
  to.cache_entry->valid = true;
 }

};

}
