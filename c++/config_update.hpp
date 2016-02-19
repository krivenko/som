#pragma once

#include <vector>
#include <utility>
#include <triqs/utility/exceptions.hpp>

#include "rectangle.hpp"
#include "configuration.hpp"

#define CHECK_CONFIG_UPDATE_OPERATIONS

namespace som {

// Object to hold a (proposed) update to a configuration
struct config_update {

 // List of changed rectangles with their indices
 // 0 <= index < INT_MAX: replace rectangle at index
 // index == INT_MAX: add rectangle
 // index < 0: remove rectanle at -index.
 std::vector<std::pair<int,rectangle>> changed_rects;

 config_update() { changed_rects.reserve(2); }

 void add_rectangle(rectangle const& r) {
  changed_rects.emplace_back(INT_MAX, r);
 }
 void remove_rectangle(int index) {
  changed_rects.emplace_back(-index-1, rectangle{});
 }
 void change_rectangle(int index, rectangle const& r) {
  changed_rects.emplace_back(index, r);
 }

 void reset() { changed_rects.resize(0); }

 friend void apply(configuration & c, config_update const& cu) {
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(cu.changed_rects.size() > 2)
   TRIQS_RUNTIME_ERROR << "config_update: only 2 operations in a row are allowed";
#endif
  for(auto const& ch : cu.changed_rects) {
   int index = ch.first;
   if(index == INT_MAX) { // add rectangle
    c.insert(ch.second);
   } else if(index < 0) { // remove rectangle
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(-index-1 >= c.size())
   TRIQS_RUNTIME_ERROR << "config_update (remove): index out of range, index = " << -index-1;
#endif
    c.remove(-index-1);
   } else { // change rectangle
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(index >= c.size())
   TRIQS_RUNTIME_ERROR << "config_update (change): index out of range, index = " << index;
#endif
    c.replace(index, ch.second);
   }
  }
 }

 int size() const { return changed_rects.size(); }
 operator bool() { return changed_rects.size() > 0; }
};

}
