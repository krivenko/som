/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2025 Igor Krivenko
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

#include <climits>
#include <utility>
#include <vector>

#include "../cache_index.hpp"
#include "../config_update.hpp"
#include "../configuration.hpp"
#include "../rectangle.hpp"

#include "observables.hpp"

namespace som {

// Integral kernel
template <observable_kind, typename Mesh> class kernel;

// Base class for all kernels (CRTP)
template <typename Derived, typename ResultType> class kernel_base {

  // Precomputed LHS values
  mutable std::vector<ResultType> lhs_cache;

  inline void check_cache_size(const cache_index& ci) const {
    if(lhs_cache.size() < ci.size())
      lhs_cache.resize(ci.size(), ResultType(first_dim(lhs_cache[0])));
  }

public:
#define DERIVED static_cast<Derived const*>(this)

  explicit kernel_base(size_t res_size) {
    lhs_cache.resize(CACHE_SIZE, ResultType(res_size));
  }

  inline ResultType const& operator()(rectangle const& rect) const {
    auto const& p = rect.cache_ptr;
    auto& ce = p.deref();
    check_cache_size(p.get_ci());
    auto& res = lhs_cache[p.id];
    // Do we have a precomputed LHS for rectangle r?
    if(!ce.valid) {
      DERIVED->apply(rect, res);
      ce.valid = true;
    }
    return res;
  }

  inline ResultType const& operator()(configuration const& c) const {
    auto const& p = c.cache_ptr;
    auto& ce = p.deref();
    check_cache_size(p.get_ci());
    auto& res = lhs_cache[p.id];
    // Do we have a precomputed LHS for configuration c?
    if(!ce.valid) {
      res() = 0;
      // cppcheck-suppress useStlAlgorithm
      for(auto const& r : c) res += operator()(r);
      ce.valid = true;
    }
    return res;
  }

  inline ResultType const& operator()(config_update const& cu) const {
    auto const& conf = cu.get_config();

    auto const& p = cu.cache_ptr;
    auto& ce = p.deref();
    check_cache_size(p.get_ci());
    auto& res = lhs_cache[p.id];
    // Do we have a precomputed LHS for config_update cu?
    if(!ce.valid) {
      // Cache entry holds the LHS for the *updated configuration*
      res() = operator()(conf);

      auto rect_it = std::begin(cu.new_rects);
      for(int index : cu.changed_indices) {
        if(index == INT_MAX) { // add rectangle
          res += operator()(*rect_it);
          ++rect_it;
        } else if(index < 0) { // remove rectangle
          res -= operator()(conf[-index - 1]);
        } else { // change rectangle
          res += operator()(*rect_it) - operator()(conf[index]);
          ++rect_it;
        }
      }
      ce.valid = true;
    }
    return res;
  }

  // Apply kernel to a configuration ignoring the cache
  inline ResultType apply_wo_caching(configuration const& c) const {
    ResultType res = ResultType::zeros({first_dim(lhs_cache[0])});
    ResultType tmp(first_dim(lhs_cache[0]));
    for(auto const& r : c) {
      DERIVED->apply(r, tmp);
      res += tmp;
    }
    return res;
  }

#undef DERIVED

  // Copy a precomputed LHS value bound to one object to a cache entry bound to
  // another object. The destination cache entry will be marked as valid
  template <typename T1, typename T2>
  void cache_copy(T1 const& from, T2 const& to) const {
    check_cache_size(from.cache_ptr.get_ci());
    lhs_cache[to.cache_ptr.id]() = lhs_cache[from.cache_ptr.id];
    to.cache_ptr.deref().valid = true;
  }

  // Swap two LHS values each bound to a given object
  // The destination cache entry (second argument) will be marked as valid
  template <typename T1, typename T2>
  void cache_swap(T1 const& from, T2 const& to) const {
    check_cache_size(from.cache_ptr.get_ci());
    std::swap(lhs_cache[to.cache_ptr.id], lhs_cache[from.cache_ptr.id]);
    to.cache_ptr.deref().valid = true;
  }
};

} // namespace som
