#pragma once

#include <vector>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <triqs/utility/exceptions.hpp>

namespace som {

struct rectangle {
 double center;
 double width;
 double height;

 rectangle(double center, double width, double height) :
  center(center), width(width), height(height) {}
 double norm() const { return width*height; }
 double operator()(double x) const {
  return (x >= center-width/2 && x <= center+width/2) ? height : 0;
 }

 // stream insertion
 friend std::ostream& operator<<(std::ostream& os, rectangle r) {
  os << "(c:" << r.center << ", w:" << r.width << ", h:" << r.height << ")";
  return os;
 }
};

class configuration {

 std::vector<rectangle> rects;

public:

 configuration() = default;
 configuration(configuration const&) = default;
 configuration(std::initializer_list<rectangle> const& l) : rects(l) {}

 std::size_t size() const { return rects.size(); }
 rectangle const& operator[](std::size_t i) { return rects[i]; }

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
  for(auto& r : rects) r.height *= alpha;
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

 // constant iterator
 using const_iterator = std::vector<rectangle>::const_iterator;
 const_iterator begin() const { return rects.begin(); }
 const_iterator cbegin() const { return rects.cbegin(); }
 const_iterator end() const { return rects.end(); }
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
