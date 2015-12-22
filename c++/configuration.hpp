#pragma once

#include <vector>
#include <algorithm>
#include <iterator>
#include <complex>
#include <ostream>
#include <triqs/utility/exceptions.hpp>
#include <triqs/arrays/vector.hpp>

namespace som {

using namespace triqs::arrays;

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

 std::complex<double> hilbert_transform(std::complex<double> z) const {
  return -height*std::log((center + width/2 - z)/(center - width/2 - z));
 }
 vector<double> tail_coefficients(long order_min, long order_max) const {
  vector<double> data(order_max - order_min + 1);

  double e1 = center - width/2, e2 = center + width/2;
  double e1n = 1, e2n = 1;
  for(long n = order_min; n <= order_max; ++n) {
   if(n < 1) { data(n-order_min) = 0; continue; }
   e1n *= e1; e2n *= e2;
   data(n-order_min) = height*(e2n - e1n)/n;
  }

  return data;
 }

 // stream insertion
 friend std::ostream& operator<<(std::ostream& os, rectangle r) {
  os << "(c:" << r.center << ", w:" << r.width << ", h:" << r.height << ")";
  return os;
 }
};

class configuration {

 std::vector<rectangle> rects;

 static const std::size_t default_max_rects = 70;

public:

 configuration(std::size_t reserved_rects = default_max_rects) {
  rects.reserve(reserved_rects);
 }
 configuration(configuration const&) = default;
 configuration(std::initializer_list<rectangle> const& l) : rects(l) {
  rects.reserve(default_max_rects);
 }

 // Number of rectangles
 std::size_t size() const { return rects.size(); }
 // Maximum number of rectangles
 std::size_t max_size() const { return rects.capacity(); }

 // Access a rectangle by index
 rectangle const& operator[](std::size_t i) const { return rects[i]; }
 rectangle & operator[](std::size_t i) { return rects[i]; }

 // Insert a new rectangle
 rectangle& insert(rectangle const& r) {
  rects.push_back(r);
  return rects.back();
 }

 // Remove a rectangle by index
 void remove(std::size_t i) {
  std::swap(rects[i],rects.back());
  rects.pop_back();
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
  for(auto & r : rects) r.height *= alpha;
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
  for(auto & r : rects) r.height *= norm / old_norm;
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
