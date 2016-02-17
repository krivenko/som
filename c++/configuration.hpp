#pragma once

#include <utility>
#include <vector>
#include <stack>
#include <algorithm>
#include <iterator>
#include <complex>
#include <ostream>
#include <triqs/utility/exceptions.hpp>
#include <triqs/arrays/vector.hpp>

#define RECT_IDS    1024
#define CONFIG_IDS  1024
#define CHECK_CONFIG_UPDATE_OPERATIONS

namespace som {

using namespace triqs::arrays;

class rectangle {
 double c;                  // center
 double w;                  // width
 double h;                  // height
 int id;                    // unique id of this rectangle

 // List of rectangle ID's free for use
 static std::stack<int> free_ids;

public:

 mutable bool lhs_computed = false; // have we already applied the kernel to this rectangle?

 static void allocate_ids(int size = RECT_IDS) {
  for(int i = size - 1; i >= 0; --i) free_ids.push(i);
 }

 rectangle() = default;
 rectangle(double center, double width, double height) :
  id(free_ids.top()), c(center), w(width), h(height) { free_ids.pop(); }
 rectangle(rectangle const& r) : id(free_ids.top()), c(r.c), w(r.w), h(r.h) { free_ids.pop(); }
 rectangle(rectangle &&) = default;
 rectangle & operator=(rectangle const&) = default;
 rectangle & operator=(rectangle &&) = default;

 ~rectangle() { free_ids.push(id); }

 double center() const { return c; }
 double width() const { return w; }
 double height() const { return h; }

  // Returns configuration id
 int get_id() const { return id; }

 double norm() const { return w*h; }
 double operator()(double x) const {
  return (x >= c-w/2 && x <= c+w/2) ? h : 0;
 }

 std::complex<double> hilbert_transform(std::complex<double> z) const {
  return -h*std::log((c + w/2 - z)/(c - w/2 - z));
 }
 vector<double> tail_coefficients(long order_min, long order_max) const {
  vector<double> data(order_max - order_min + 1);

  double e1 = c - w/2, e2 = c + w/2;
  double e1n = 1, e2n = 1;
  for(long n = order_min; n <= order_max; ++n) {
   if(n < 1) { data(n-order_min) = 0; continue; }
   e1n *= e1; e2n *= e2;
   data(n-order_min) = h*(e2n - e1n)/n;
  }

  return data;
 }

 // stream insertion
 friend std::ostream& operator<<(std::ostream& os, rectangle r) {
  os << "(c:" << r.c << ", w:" << r.w << ", h:" << r.h << ")";
  return os;
 }

 rectangle const& operator*=(double alpha) { h *= alpha; return *this; }
};
std::stack<int> rectangle::free_ids;

class configuration {

 // Rectangles in this configuration
 std::vector<rectangle> rects;
 // Unique id of this configuration
 int id;

 // List of rectangle ID's free for use
 static std::stack<int> free_ids;

 static const int default_max_rects = 70;

public:

 mutable bool lhs_computed = false; // have we already applied the kernel to this configuration?

 static void allocate_ids(int size = CONFIG_IDS) {
  for(int i = size - 1; i >= 0; --i) free_ids.push(i);
 }

 configuration(int reserved_rects = default_max_rects) : id(free_ids.top()) {
  free_ids.pop();
  rects.reserve(reserved_rects);
 }
 configuration(std::initializer_list<rectangle> const& l) : id(free_ids.top()), rects(l) {
  free_ids.pop();
  rects.reserve(default_max_rects);
 }
 configuration(configuration const& r) : id(free_ids.top()), rects(r.rects) { free_ids.pop(); }
 configuration(configuration &&) = default;
 configuration & operator=(configuration const&) = default;
 configuration & operator=(configuration &&) = default;

 ~configuration() { free_ids.push(id); }

 // Number of rectangles
 int size() const { return rects.size(); }
 // Maximum number of rectangles
 int max_size() const { return rects.capacity(); }
 // Returns configuration id
 int get_id() const { return id; }

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
  for(auto & r : rects) r *= alpha;
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
  for(auto & r : rects) r *= norm / old_norm;
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
std::stack<int> configuration::free_ids;

// Object to hold a (proposed) update to a configuration
class config_update {

 // List of changed rectangles with their indices
 // 0 <= index < INT_MAX: replace rectangle at index
 // index == INT_MAX: add rectangle
 // index < 0: remove rectanle at -index.
 std::vector<std::pair<int,rectangle>> changed_rects;

public:

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

 void apply(configuration & c) {
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(changed_rects.size() > 2)
   TRIQS_RUNTIME_ERROR << "config_update: only 2 operations in a row are allowed";
#endif
  for(auto const& ch : changed_rects) {
   int index = ch.first;
   if(index == INT_MAX) { // add rectangle
    c.insert(ch.second);
   } else if(index < 0) { // remove rectangle
#ifdef CHECK_CONFIG_UPDATE_OPERATIONS
  if(-index >= c.size())
   TRIQS_RUNTIME_ERROR << "config_update (remove): index out of range, index = " << -index;
#endif
    c.remove(-index);
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
