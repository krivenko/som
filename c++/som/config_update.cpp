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

#include <climits>
#include <utility>

#include <triqs/utility/exceptions.hpp>

#include "config_update.hpp"

#define CHECK_CONFIG_UPDATE_OPERATIONS

namespace som {

config_update::config_update(configuration& conf, cache_index& ci)
   : conf(conf), ci(ci), cache_id(ci.acquire()) {
  changed_indices.reserve(2);
  new_rects.reserve(2);
}
config_update::config_update(config_update const& cu)
   : conf(cu.conf)
   , changed_indices(cu.changed_indices)
   , new_rects(cu.new_rects)
   , ci(cu.ci)
   , cache_id(ci.acquire()) {}
config_update::config_update(config_update&& cu) noexcept
   : conf(cu.conf)
   , changed_indices(cu.changed_indices)
   , new_rects(cu.new_rects)
   , ci(cu.ci)
   , cache_id(ci.acquire()) {}
config_update::~config_update() { ci.decref(cache_id); }

void config_update::add_rectangle(rectangle const& r) {
  changed_indices.push_back(INT_MAX);
  new_rects.push_back(r);
  ci[cache_id].valid = false;
}
void config_update::add_rectangle(rectangle&& r) {
  changed_indices.push_back(INT_MAX);
  new_rects.push_back(std::move(r));
  ci[cache_id].valid = false;
}

void config_update::remove_rectangle(int index) {
  changed_indices.push_back(-index - 1);
  ci[cache_id].valid = false;
}

void config_update::change_rectangle(int index, rectangle const& r) {
  changed_indices.push_back(index);
  new_rects.push_back(r);
  ci[cache_id].valid = false;
}
void config_update::change_rectangle(int index, rectangle&& r) {
  changed_indices.push_back(index);
  new_rects.push_back(std::move(r));
  ci[cache_id].valid = false;
}

void config_update::reset() {
  changed_indices.clear();
  new_rects.clear();
  ci[cache_id].valid = false;
}

void config_update::apply() {
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(changed_indices.size() > 2)
    TRIQS_RUNTIME_ERROR
        << "config_update: only 2 operations in a row are allowed";
#endif
  auto rect_it = std::begin(new_rects);
  for(auto index : changed_indices) {
    if(index == INT_MAX) { // add rectangle
      conf.insert(*rect_it);
      ++rect_it;
    } else if(index < 0) { // remove rectangle
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
      if(-index - 1 >= conf.size())
        TRIQS_RUNTIME_ERROR
            << "config_update (remove): index out of range, index = "
            << -index - 1;
#endif
      conf.remove(-index - 1);
    } else { // change rectangle
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
      if(index >= conf.size())
        TRIQS_RUNTIME_ERROR
            << "config_update (change): index out of range, index = " << index;
#endif
      conf.replace(index, *rect_it);
      ++rect_it;
    }
  }
  ci[conf.cache_id].valid = false;
  reset();
}

} // namespace som
