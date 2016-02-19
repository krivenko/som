#pragma once

#include <vector>
#include <complex>
#include <ostream>
#include <triqs/arrays/vector.hpp>

#include "ids.hpp"

namespace som {

using namespace triqs::arrays;

class rectangle {
 double c;                  // center
 double w;                  // width
 double h;                  // height
 int id;                    // unique id of this rectangle

public:

 mutable bool lhs_computed = false; // have we already applied the kernel to this rectangle?

 rectangle() = default;
 rectangle(double center, double width, double height) :
  id(get_rectangle_id()), c(center), w(width), h(height) {}
 rectangle(rectangle const& r) : id(get_rectangle_id()), c(r.c), w(r.w), h(r.h) {}
 rectangle(rectangle &&) = default;
 rectangle & operator=(rectangle const&) = default;
 rectangle & operator=(rectangle &&) = default;

 ~rectangle() { release_rectangle_id(id); }

 inline double center() const { return c; }
 inline double width() const { return w; }
 inline double height() const { return h; }

  // Returns configuration id
 inline int get_id() const { return id; }

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

}
