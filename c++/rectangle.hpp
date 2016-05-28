/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016 by I. Krivenko
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

 // Pointer to a cache entry descriptor
 // Rectangles are immutable objects, therefore they acquire a new cache entry
 // descriptor only on construction. Copies inherit the same descriptor, and
 // increase the refcount to it. The refcount is decreased by one on destruction.
 // Cache entry bound to an existing rectangle is never invalidated, once the LHS is computed.
 cache_index::entry * cache_entry;

public:

 rectangle(double center, double width, double height, cache_index & ci) :
  center(center), width(width), height(height), cache_entry(ci.aquire()) {}

 rectangle(rectangle const& r) :
  center(r.center), width(r.width), height(r.height), cache_entry(r.cache_entry) {
  cache_entry->incref();
 }
 rectangle(rectangle && r) :
  center(r.center), width(r.width), height(r.height), cache_entry(r.cache_entry) {
  cache_entry->incref();
 }
 rectangle & operator=(rectangle const& r) {
  cache_entry->decref();
  center = r.center; width = r.width; height = r.height; cache_entry = r.cache_entry;
  cache_entry->incref();
  return *this;
 }
 rectangle & operator=(rectangle && r) {
  using std::swap;
  swap(center, r.center);
  swap(width, r.width);
  swap(height, r.height);
  swap(cache_entry, r.cache_entry);
  return *this;
 }
 ~rectangle() { cache_entry->decref(); }

 bool operator==(rectangle const& r) const {
  using triqs::utility::is_zero;
  return is_zero(center-r.center, 1e-9) &&
         is_zero(width-r.width, 1e-9) &&
         is_zero(height-r.height, 1e-9);
 }
 bool operator!=(rectangle const& r) const { return !operator==(r); }

 double norm() const { return width * height; }
 double operator()(double x) const {
  return (x >= center - width/2 && x <= center + width/2) ? height : 0;
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
 friend std::ostream& operator<<(std::ostream& os, rectangle const& r) {
  os << "(c:" << r.center << ", w:" << r.width << ", h:" << r.height << ")";
  return os;
 }

 // multiplication by scalar
 friend rectangle operator*(rectangle const& r, double alpha) {
  return {r.center, r.width, r.height * alpha, r.cache_entry->index};
 }
 friend rectangle operator*(double alpha, rectangle const& r) { return r*alpha; }

 // POD version of the rectangle, used in MPI operations
 struct pod_t {
  double center; double width; double height;
  pod_t() = default;
  pod_t(rectangle const& r) : center(r.center), width(r.width), height(r.height) {}
 };

};

}

namespace triqs { namespace mpi {

template<> inline MPI_Datatype mpi_datatype<som::rectangle::pod_t>() {
 bool type_committed = false;
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
