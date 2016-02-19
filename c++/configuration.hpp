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

#include <utility>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <triqs/utility/exceptions.hpp>

#include "rectangle.hpp"

namespace som {

using namespace triqs::arrays;

class configuration {

 // Rectangles in this configuration
 std::vector<rectangle> rects;

 static const int default_max_rects = 70;

public:

 configuration(int reserved_rects = default_max_rects) {
  rects.reserve(reserved_rects);
 }
 configuration(std::initializer_list<rectangle> const& l) : rects(l) {
  rects.reserve(default_max_rects);
 }
 configuration(configuration const&) = default;
 //configuration(configuration &&) = default;
 configuration & operator=(configuration const&) = default;
 //configuration & operator=(configuration &&) = default;

 // Number of rectangles
 int size() const { return rects.size(); }
 // Maximum number of rectangles
 int max_size() const { return rects.capacity(); }
 // Norm of this configuration
 double norm() const {
  double res = 0;
  for(auto const& r : rects) res += r.norm();
  return res;
 }

 // Access a rectangle by index
 rectangle const& operator[](int index) const { return rects[index]; }

 // Insert a new rectangle
 rectangle& insert(rectangle const& r) {
  rects.push_back(r);
  return rects.back();
 }

 // Remove a rectangle by index
 void remove(int index) {
  std::swap(rects[index],rects.back());
  rects.pop_back();
 }

 // Replace a rectangle
 void replace(int index, rectangle const& r) {
  rects[index] = r;
 }

 // Sum of configurations: all rectangles from both of them
 configuration& operator+=(configuration const& c) {
  rects.reserve(rects.size() + c.rects.size());
  std::copy(c.rects.begin(),c.rects.end(),std::back_inserter(rects));
  return *this;
 }
 configuration operator+(configuration const& c) const {
  configuration s(*this);
  s += c;
  return s;
 }
 // Multiply configuration by a positive scalar
 // Heights of all rectangles are scaled
 configuration& operator*=(double alpha) {
  if(alpha < 0) TRIQS_RUNTIME_ERROR <<
                "Cannot multiply a configuration by a negative number " << alpha;
  for(auto & r : rects) r = r*alpha;
  return *this;
 }
 friend configuration operator*(configuration const& c, double alpha) {
  configuration p(c);
  p *= alpha;
  return p;
 }
 friend configuration operator*(double alpha, configuration const& c) {
  return c*alpha;
 }

 // Normalize configuration to have a total area of norm
 void normalize(double norm = 1.0) {
  if(norm <= 0) TRIQS_RUNTIME_ERROR <<
                "Configuration must have a positive norm; norm = "
                << norm << " is requested";
  double old_norm = 0;
  for(auto const& r : rects) old_norm += r.norm();
  for(auto & r : rects) r = r * (norm / old_norm);
 }

 // constant iterator
 using const_iterator = std::vector<rectangle>::const_iterator;
 const_iterator begin() const { return rects.begin(); }
 const_iterator end() const { return rects.end(); }
 const_iterator cbegin() const { return rects.cbegin(); }
 const_iterator cend() const { return rects.cend(); }

 // stream insertion
 friend std::ostream& operator<<(std::ostream& os, configuration c) {
  bool print_comma = false;
  for(auto const& r : c.rects) {
   if(print_comma) os << ",";
   else print_comma = true;
   os << r;
  }
  return os;
 }
};

}
