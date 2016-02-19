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

#define CHECK_CONFIG_UPDATE_OPERATIONS

namespace som {

// Object to hold a (proposed) update to a configuration
struct config_update {

 // configuration this update will be applied to
 configuration & conf;

 // List of changed rectangles with their indices
 // 0 <= index < INT_MAX: replace rectangle at index
 // index == INT_MAX: add rectangle
 // index < 0: remove rectanle at -index.
 std::vector<std::pair<int,rectangle>> changed_rects;

 config_update(configuration & conf) : conf(conf) { changed_rects.reserve(2); }

 void add_rectangle(rectangle const& r) {
  changed_rects.emplace_back(INT_MAX, r);
 }
 void remove_rectangle(int index) {
  changed_rects.emplace_back(-index-1, rectangle{});
 }
 void change_rectangle(int index, rectangle const& r) {
  changed_rects.emplace_back(index, r);
 }

 // access the base configuration
 configuration const& get_config() const { return conf; }

 void reset() { changed_rects.resize(0); }

 void apply() {
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(changed_rects.size() > 2)
   TRIQS_RUNTIME_ERROR << "config_update: only 2 operations in a row are allowed";
#endif
  for(auto const& ch : changed_rects) {
   int index = ch.first;
   if(index == INT_MAX) { // add rectangle
    conf.insert(ch.second);
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
    conf.replace(index, ch.second);
   }
  }
  reset();
 }

 int size() const { return changed_rects.size(); }
 operator bool() { return changed_rects.size() > 0; }
};

}
