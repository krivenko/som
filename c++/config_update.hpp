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
#include <utility>
#include <triqs/utility/exceptions.hpp>

#include "configuration.hpp"
#include "cache_index.hpp"

#define CHECK_CONFIG_UPDATE_OPERATIONS

namespace som {

// Object to hold a (proposed) update to a configuration
struct config_update {

 // configuration this update will be applied to
 configuration & conf;

 // List of indices of changed rectangles
 // 0 <= index < INT_MAX: replace rectangle at index
 // index == INT_MAX: add rectangle
 // index < 0: remove rectanle at -index.
 std::vector<int> changed_indices;
 // New rectangles (for insertions and chages)
 std::vector<rectangle> new_rects;

 // Pointer to a cache entry descriptor
 // Every new config_update object, including copies, aquire a new
 // cache entry descriptor, which is then released in the destructor.
 // Methods reset() and *_rectangle() will invalidate the cache entry.
 cache_index::entry * cache_entry;

 config_update(configuration & conf, cache_index & ci) :
  conf(conf), cache_entry(ci.aquire()) {
  changed_indices.reserve(2);
  new_rects.reserve(2);
 }
 config_update(config_update const& cu) :
  conf(cu.conf), changed_indices(cu.changed_indices), new_rects(cu.new_rects),
  cache_entry(cu.cache_entry->index.aquire()) {}
 config_update(config_update && cu) :
  conf(cu.conf), changed_indices(std::move(cu.changed_indices)), new_rects(std::move(cu.new_rects)),
  cache_entry(cu.cache_entry->index.aquire()) {}
 config_update & operator=(config_update const&) = delete;
 config_update & operator=(config_update &&) = delete;
 ~config_update() { cache_entry->decref(); }

 void add_rectangle(rectangle const& r) {
  changed_indices.push_back(INT_MAX);
  new_rects.push_back(r);
  cache_entry->valid = false;
 }
 void add_rectangle(rectangle && r) {
  changed_indices.push_back(INT_MAX);
  new_rects.push_back(std::move(r));
  cache_entry->valid = false;
 }

 void remove_rectangle(int index) {
  changed_indices.push_back(-index-1);
  cache_entry->valid = false;
 }

 void change_rectangle(int index, rectangle const& r) {
  changed_indices.push_back(index);
  new_rects.push_back(r);
  cache_entry->valid = false;
 }
 void change_rectangle(int index, rectangle && r) {
  changed_indices.push_back(index);
  new_rects.push_back(std::move(r));
  cache_entry->valid = false;
 }

 // access the base configuration
 configuration const& get_config() const { return conf; }

 void reset() {
  changed_indices.clear();
  new_rects.clear();
  cache_entry->valid = false;
 }

 void apply() {
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(changed_indices.size() > 2)
   TRIQS_RUNTIME_ERROR << "config_update: only 2 operations in a row are allowed";
#endif
  auto rect_it = std::cbegin(new_rects);
  for(auto index : changed_indices) {
   if(index == INT_MAX) { // add rectangle
    conf.insert(*rect_it);
    ++rect_it;
   } else if(index < 0) { // remove rectangle
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(-index-1 >= conf.size())
   TRIQS_RUNTIME_ERROR << "config_update (remove): index out of range, index = " << -index-1;
#endif
    conf.remove(-index-1);
   } else { // change rectangle
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(index >= conf.size())
   TRIQS_RUNTIME_ERROR << "config_update (change): index out of range, index = " << index;
#endif
    conf.replace(index, *rect_it);
    ++rect_it;
   }
  }
  conf.cache_entry->valid = false;
  reset();
 }

 int size() const { return changed_indices.size(); }
 operator bool() { return changed_indices.size() > 0; }
};

}
