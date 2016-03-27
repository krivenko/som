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

#include <vector>
#include <stack>
#include <utility>

namespace som {

// cache_index establishes connection between rectangles/configurations/configuration updates
// and the internal LHS cache stored inside kernel_base.
//
// See rectangle.hpp, configuration.hpp and config_update.hpp for details
// on how the respective classes interact with this object
class cache_index {

public:

 // Descriptor of a cache entry
 struct entry {
  const int id;     // position in the cache
  int refcount;     // number of references to the respective cache entry
  bool valid;       // is the respective cache entry valid?

  // Reference to the parent cache_index object
  cache_index & index;

  entry(int id, cache_index & index) : id(id), refcount(0), valid(false), index(index) {}

  // Increase reference count
  inline void incref() { ++refcount; }
  // Decrease reference count
  inline void decref() { --refcount; if(refcount == 0) index.spare_ids.push(id); }
 };

private:

 std::vector<entry> entries; // Cache entry descriptors
 std::stack<int> spare_ids;  // Pool of descriptors not referenced by any objects

public:

 // Constructor
 cache_index(int capacity = CACHE_SIZE) {
  entries.reserve(capacity);
  for(int id = 0; id < capacity; ++id) {
   entries.emplace_back(id, *this);
   spare_ids.push(capacity - 1 - id);
  }
 }

 // Aquire ownership over a free cache entry
 entry * aquire() {
  if(spare_ids.size() == 0)
   TRIQS_RUNTIME_ERROR << "Capacity of the LHS cache is exhausted. "
                          "Consider rebuilding with larger CACHE_SIZE.";
  int id = spare_ids.top();
  spare_ids.pop();
  auto h = &entries[id];
  ++h->refcount;
  h->valid = false;
  return h;
 }

 // Mark all cache entries as invalid
 void invalidate_all() { for(auto & h : entries) h.valid = false; }

};

}
