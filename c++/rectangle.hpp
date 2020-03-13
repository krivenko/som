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
#include <complex>
#include <ostream>

#include <triqs/arrays/vector.hpp>
#include <triqs/utility/numeric_ops.hpp>
#include <triqs/mpi/base.hpp>

#include "cache_index.hpp"

namespace som {

using namespace triqs::arrays;

struct rectangle {

 double center;             // center
 double width;              // width
 double height;             // height

 // Reference to the cache index and id within the cache
 // Rectangles are immutable objects, therefore they acquire a new cache entry
 // descriptor only on construction. Copies inherit the same descriptor, and
 // increase the refcount to it. The refcount is decreased by one on destruction.
 // Cache entry bound to an existing rectangle is never invalidated, once the LHS is computed.
 cache_index & ci;
 int cache_id;

public:

 rectangle(double center, double width, double height, cache_index & ci) :
  center(center), width(width), height(height), ci(ci), cache_id(ci.aquire()) {}

 rectangle(rectangle const& r) :
  center(r.center), width(r.width), height(r.height), ci(r.ci), cache_id(r.cache_id) {
  ci.incref(cache_id);
 }
 rectangle(rectangle && r) noexcept :
  center(r.center), width(r.width), height(r.height), ci(r.ci), cache_id(r.cache_id) {
  ci.incref(cache_id);
 }
 rectangle & operator=(rectangle const& r) {
  ci.decref(cache_id);
  cache_id = r.cache_id;
  center = r.center; width = r.width; height = r.height;
  ci.incref(cache_id);
  return *this;
 }
 rectangle & operator=(rectangle && r) noexcept {
  using std::swap;
  swap(center, r.center);
  swap(width, r.width);
  swap(height, r.height);
  swap(cache_id, r.cache_id);
  return *this;
 }
 ~rectangle() { ci.decref(cache_id); }

 bool operator==(rectangle const& r) const {
  using triqs::utility::is_zero;
  return is_zero(center-r.center, 1e-8) &&
         is_zero(width-r.width, 1e-8) &&
         is_zero(height-r.height, 1e-8);
 }
 bool operator!=(rectangle const& r) const { return !operator==(r); }

 double norm() const { return width * height; }
 double operator()(double x) const {
  return (x >= center - width/2 && x <= center + width/2) ? height : 0;
 }

 std::complex<double> hilbert_transform(std::complex<double> z, bool multiply_by_e = false) const {
  if(multiply_by_e)
   // -h \int_{c-w/2}^{c+w/2} d\epsilon' \frac{\epsilon'}{\epsilon' - \epsilon - i0}
   return -height*(width + z*std::log((center + width/2 - z)/(center - width/2 - z)));
  else
   // -h \int_{c-w/2}^{c+w/2} d\epsilon' \frac{1}{\epsilon' - \epsilon - i0}
   return -height*std::log((center + width/2 - z)/(center - width/2 - z));
 }
 vector<double> tail_coefficients(long order_min, long order_max, bool multiply_by_e = false) const {
  vector<double> data(order_max - order_min + 1);

  double e1 = center - width/2, e2 = center + width/2;
  double e1n, e2n;
  int denom_shift;
  if(multiply_by_e) { e1n = e1; e2n = e2; denom_shift = 1; }
  else              { e1n = 1.0; e2n = 1.0; denom_shift = 0; }
  for(long n = order_min; n <= order_max; ++n) {
   if(n < 1) { data(n-order_min) = 0; continue; }
   e1n *= e1; e2n *= e2;
   data(n-order_min) = height*(e2n - e1n)/(n + denom_shift);
  }

  return data;
 }

 // stream insertion
 friend std::ostream& operator<<(std::ostream& os, rectangle const& r) {
  os << "(c:" << r.center << ", w:" << r.width << ", h:" << r.height << ")";
  return os;
 }

 // multiplication by scalar
 friend rectangle operator*(rectangle const& r, double alpha) {
  return {r.center, r.width, r.height * alpha, r.ci};
 }
 friend rectangle operator*(double alpha, rectangle const& r) { return r*alpha; }

 // POD version of the rectangle, used in MPI operations
 struct pod_t {
  double center; double width; double height;
  pod_t() = default;
  pod_t(rectangle const& r) : center(r.center), width(r.width), height(r.height) {}
 };

 // Convert to tuple (center,width,height)
 operator std::tuple<double,double,double>() const {
  return std::make_tuple(center,width,height);
 }

};

}

namespace triqs { namespace mpi {

template<> inline MPI_Datatype mpi_datatype<som::rectangle::pod_t>() {
 static bool type_committed = false;
 static MPI_Datatype dt;
 if(!type_committed) {
  int blocklengths[] = {1,1,1};
  MPI_Aint displacements[] = {0,sizeof(double),2*sizeof(double)};
  MPI_Datatype types[] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  MPI_Type_create_struct(3, blocklengths, displacements, types, &dt);
  MPI_Type_commit(&dt);
  type_committed = true;
 }
 return dt;
}

}}
