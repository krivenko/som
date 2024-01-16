/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <cassert>
#include <stack>
#include <utility>
#include <vector>

namespace som {

// cache_index establishes connection between
// rectangles/configurations/configuration updates and the internal LHS cache
// stored inside kernel_base.
//
// See rectangle.hpp, configuration.hpp and config_update.hpp for details
// on how the respective classes interact with this object
class cache_index {

public:
  // Descriptor of a cache entry
  struct entry {
    int refcount = 0;   // number of references to the respective cache entry
    bool valid = false; // is the respective cache entry valid?
  };

private:
  // Cache entry descriptors
  std::vector<entry> entries;
  // Pool of descriptors not referenced by any objects
  std::stack<std::size_t> spare_ids;

  void extend();

public:
  // Constructor
  cache_index() { extend(); }

  // Acquire ownership over a free cache entry
  std::size_t acquire();

  // Increase reference count for cache entry id
  inline void incref(std::size_t id) { ++entries[id].refcount; }

  // Decrease reference count
  inline void decref(std::size_t id) {
    auto& e = entries[id];
    --e.refcount;
    if(e.refcount == 0) spare_ids.push(id);
  }

  // Access to entries
  inline entry& operator[](std::size_t id) { return entries[id]; }
  inline entry const& operator[](std::size_t id) const { return entries[id]; }

  // Current size of the index
  [[nodiscard]] inline std::size_t size() const { return entries.size(); }

  // Mark all cache entries as invalid
  void invalidate_all();

  // Number of currently used cache entries
  [[nodiscard]] std::size_t n_used_entries() const {
    return entries.size() - spare_ids.size();
  }
};

// Pointer to a cache entry
struct cache_entry_ptr {
  cache_index* ci = nullptr; // Cache index the entry is part of
  std::size_t id = 0;        // ID of the entry

  // Construct a null pointer
  cache_entry_ptr() = default;
  // Construct a valid pointer
  inline explicit cache_entry_ptr(cache_index& ci)
     : ci(&ci), id(ci.acquire()) {}
  // Destroy the pointer
  inline ~cache_entry_ptr() {
    if(ci) ci->decref(id);
  }
  // Copy-constructor
  inline cache_entry_ptr(cache_entry_ptr const& p) : ci(p.ci), id(p.id) {
    if(ci) ci->incref(id);
  }
  // Move-constructor
  inline cache_entry_ptr(cache_entry_ptr&& p) noexcept : ci(p.ci), id(p.id) {
    if(ci) ci->incref(id);
  }
  // Copy-assignment
  inline cache_entry_ptr& operator=(cache_entry_ptr const& p) {
    if(ci) ci->decref(id);
    ci = p.ci;
    id = p.id;
    if(ci) ci->incref(id);
    return *this;
  }
  // Move-assignment
  inline cache_entry_ptr& operator=(cache_entry_ptr&& p) noexcept {
    using std::swap;
    swap(ci, p.ci);
    swap(id, p.id);
    return *this;
  }

  // Null pointer?
  inline operator bool() const { return ci != nullptr; }

  // Dereference this pointer
  [[nodiscard]] inline cache_index::entry& deref() const {
    assert(ci);
    return (*ci)[id];
  }

  // Dereference 'ci'
  [[nodiscard]] inline cache_index& get_ci() const {
    assert(ci);
    return *ci;
  }

  // Mark the pointee as invalid
  inline void invalidate_entry() {
    if(ci) (*ci)[id].valid = false;
  }
};

} // namespace som
