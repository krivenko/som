#pragma once

#include <vector>

#include "../configuration.hpp"

namespace som {

// Kinds of observables
enum observable_kind {FermionicGf, Susceptibility, Conductivity};

// Integral kernels
template<observable_kind Kind, typename Mesh> class kernel;

using namespace triqs::arrays;

// Container for precomputed (kernel*rectangle) and (kernel*configuration) vectors
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

 ResultType const& operator()(rectangle const& rect) const {

  ResultType & res = rect_cache[rect.get_id()];
  if(rect.lhs_computed) return res; // Satisfied from cache
  static_cast<Derived const*>(this)->apply(rect, res);
  rect.lhs_computed = true;
  return res;
 }

 ResultType const& operator()(configuration const& c) const {

  ResultType & res = config_cache[c.get_id()];
  if(c.lhs_computed) return res;  // Satisfied from cache
  for(auto const& r : c) res += *static_cast<Derived*>(this)(r);
  c.lhs_computed = true;
  return res;
 }

 // Multiply cached values by a constant
 void cache_multiply(rectangle const& rect, double alpha) {
  rect_cache[rect.get_id()] *= alpha;
 }
 void cache_multiply(configuration const& config, double alpha) {
  config_cache[config.get_id()] *= alpha;
 }
};

}
