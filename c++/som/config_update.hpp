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

#pragma once

#include <vector>

#include "cache_index.hpp"
#include "configuration.hpp"
#include "rectangle.hpp"

#define CHECK_CONFIG_UPDATE_OPERATIONS

namespace som {

// Object to hold a (proposed) update to a configuration
struct config_update {

  // configuration this update will be applied to
  configuration& conf;

  // List of indices of changed rectangles
  // 0 <= index < INT_MAX: replace rectangle at index
  // index == INT_MAX: add rectangle
  // index < 0: remove rectanle at -index.
  std::vector<int> changed_indices;
  // New rectangles (for insertions and chages)
  std::vector<rectangle> new_rects;

  // Reference to the cache index and id within the cache
  // Every new config_update object, including copies, acquire a new
  // cache entry descriptor, which is then released in the destructor.
  // Methods reset() and *_rectangle() will invalidate the cache entry.
  cache_index& ci;
  int cache_id;

  config_update(configuration& conf, cache_index& ci);
  config_update(config_update const& cu);
  config_update(config_update&& cu) noexcept;
  config_update& operator=(config_update const&) = delete;
  config_update& operator=(config_update&&) = delete;
  ~config_update();

  void add_rectangle(rectangle const& r);
  void add_rectangle(rectangle&& r);

  void remove_rectangle(int index);

  void change_rectangle(int index, rectangle const& r);
  void change_rectangle(int index, rectangle&& r);

  // Access the base configuration
  [[nodiscard]] configuration const& get_config() const { return conf; }

  void reset();

  void apply();

  [[nodiscard]] size_t size() const { return changed_indices.size(); }
  operator bool() { return changed_indices.size() > 0; }
};

} // namespace som
