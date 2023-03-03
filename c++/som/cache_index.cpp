/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#ifdef EXT_DEBUG
#include <iostream>
#endif

#include "cache_index.hpp"

namespace som {

void cache_index::extend() {
  std::size_t cap = entries.size();
  std::size_t new_cap = cap + CACHE_SIZE;
#ifdef EXT_DEBUG
  std::cerr << "Extending LHS cache from " << cap << " to " << new_cap
            << " entries." << std::endl;
#endif
  entries.reserve(new_cap);
  for(std::size_t id = cap; id < new_cap; ++id) {
    entries.emplace_back();
    spare_ids.push(new_cap - 1 + cap - id);
  }
}

// Acquire ownership over a free cache entry
std::size_t cache_index::acquire() {
  if(spare_ids.empty()) extend();
  std::size_t id = spare_ids.top();
  spare_ids.pop();
  auto& e = entries[id];
  ++e.refcount;
  e.valid = false;
  return id;
}

void cache_index::invalidate_all() {
  for(auto& h : entries) h.valid = false;
}

} // namespace som
