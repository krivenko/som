#pragma once

#include <vector>

#include "../rectangle.hpp"
#include "../configuration.hpp"
#include "../config_update.hpp"

namespace som {

// Kinds of observables
enum observable_kind {FermionGf, Susceptibility, Conductivity};

// Integral kernels
template<observable_kind, typename Mesh> class kernel;

// Base class for all kernels (CRTP)
template<typename Derived, typename ResultType>
class kernel_base {

protected:

 // Container for precomputed vectors (kernel*rectangle)
 mutable std::vector<ResultType> rect_cache;
 // Container for precomputed vectors (kernel*configuration)
 mutable std::vector<ResultType> config_cache;

public:

 kernel_base(ResultType const& init_lhs,
             int rect_cache_size = RECT_IDS,
             int config_cache_size = CONFIG_IDS) :
  rect_cache(rect_cache_size, init_lhs), config_cache(config_cache_size, init_lhs) {}

 inline ResultType const& operator()(rectangle const& rect) const {

  ResultType & res = rect_cache[rect.get_id()];
  if(rect.lhs_computed) return res; // Satisfied from cache
  static_cast<Derived const*>(this)->apply(rect, res);
  rect.lhs_computed = true;
  return res;
 }

 inline ResultType const& operator()(configuration const& c) const {

  ResultType & res = config_cache[c.get_id()];
  if(c.lhs_computed) return res;  // Satisfied from cache
  res() = 0;
  for(auto const& r : c) res += operator()(r);
  c.lhs_computed = true;
  return res;
 }

 inline ResultType operator()(configuration const& c, config_update const& cu) const {
  ResultType res = operator()(c);
  for(auto const& ch : cu.changed_rects) {
   int index = ch.first;
   if(index == INT_MAX) { // add rectangle
    res += operator()(ch.second);
   } else if(index < 0) { // remove rectangle
    res -= operator()(c[-index-1]);
   } else { // change rectangle
    res += operator()(ch.second) - operator()(c[index]);
   }
  }
  return res;
 }

 // Multiply cached values by a constant
 inline void cache_multiply(rectangle const& rect, double alpha) {
  if(rect.lhs_computed) rect_cache[rect.get_id()] *= alpha;
 }
 inline void cache_multiply(configuration const& c, double alpha) {
  for(auto const& r : c) cache_multiply(r, alpha);
  if(c.lhs_computed) config_cache[c.get_id()] *= alpha;
 }

 // Copy cached values
 inline void cache_copy(rectangle const& from, rectangle const& to) {
  if(from.lhs_computed) {
   rect_cache[to.get_id()] = rect_cache[from.get_id()];
   to.lhs_computed = true;
  }
 }
 inline void cache_copy(configuration const& from, configuration const& to) {
  if(from.lhs_computed) {
   for(int i = 0; i < from.size(); ++i) cache_copy(from[i], to[i]);
   config_cache[to.get_id()] = config_cache[from.get_id()];
   to.lhs_computed = true;
  }
 }

 // Update cached values
 void update_config_cache(configuration const& c, config_update const& cu) const {
  config_cache[c.get_id()] = operator()(c, cu);
  c.lhs_computed = true;
 }
};

}
