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
#include <triqs/arrays.hpp>
#include <triqs/mpi/vector.hpp>

#include "rectangle.hpp"
#include "cache_index.hpp"

namespace som {

using namespace triqs::arrays;
namespace h5 = triqs::h5;

struct configuration {

 // Rectangles in this configuration
 std::vector<rectangle> rects;

 static const int default_max_rects = 70;

 // Pointer to a cache entry descriptor
 // Every new configuration object, including copies, aquire a new
 // cache entry descriptor, which is then released in the destructor.
 // Assignments and compound operations *=, += will invalidate the
 // cache entry. One should manually call kernel_base::cache_* methods
 // to update the invalidated entries without a full recomputation
 cache_index::entry * cache_entry;

private:

 // Configurations are supposed to be modified
 // only by config_update and solution_worker
 friend class config_update;
 template<typename KernelType> friend class solution_worker;

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

 // Remove all rectangles, only for internal use in h5_read
 void clear() { rects.clear(); }

public:

 configuration(cache_index & ci, int reserved_rects = default_max_rects) :
  cache_entry(ci.aquire()) {
  rects.reserve(reserved_rects);
 }
 configuration(std::initializer_list<rectangle> const& l, cache_index & ci) :
  rects(l), cache_entry(ci.aquire()) {
  rects.reserve(default_max_rects);
 }
 configuration(configuration const& c) :
  rects(c.rects), cache_entry(c.cache_entry->index.aquire()) {}
 configuration(configuration && c) :
  rects(std::move(c.rects)), cache_entry(c.cache_entry->index.aquire()) {}
 configuration & operator=(configuration const& c) {
  rects = c.rects;
  cache_entry->valid = false;
  return *this;
 }
 configuration & operator=(configuration && c) {
  std::swap(rects, c.rects);
  cache_entry->valid = false;
  return *this;
 }
 ~configuration() { cache_entry->decref(); }

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

 // Equality
 bool operator==(configuration const& c) const { return rects == c.rects; }
 bool operator!=(configuration const& c) const { return rects != c.rects; }

 // Sum of configurations: all rectangles from both of them
 configuration& operator+=(configuration const& c) {
  rects.reserve(rects.size() + c.rects.size());
  std::copy(c.rects.begin(),c.rects.end(),std::back_inserter(rects));
  cache_entry->valid = false;
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
  cache_entry->valid = false;
  return *this;
 }
 friend configuration operator*(configuration const& c, double alpha) {
  configuration p(c);
  p *= alpha;
  return p;
 }
 friend configuration operator*(double alpha, configuration const& c) {
  configuration p(c);
  p *= alpha;
  return p;
 }

 // Normalize configuration to have a total area of norm
 void normalize(double norm = 1.0) {
  if(norm <= 0) TRIQS_RUNTIME_ERROR <<
                "Configuration must have a positive norm; norm = "
                << norm << " is requested";
  double old_norm = 0;
  for(auto const& r : rects) old_norm += r.norm();
  for(auto & r : rects) r = r * (norm / old_norm);
  // Invalidate LHS cache entry
  cache_entry->valid = false;
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

 // MPI reduce
 friend configuration mpi_reduce(configuration const& c,
                                 triqs::mpi::communicator comm = {}, int root = 0, bool all = false) {
  if(comm.size() == 1) return c;
  configuration res(c.cache_entry->index);

  std::vector<rectangle::pod_t> pod_rects(std::begin(c), std::end(c));
  pod_rects = mpi_gather(pod_rects, comm, root, all, std::true_type());

  for(auto const& r : pod_rects) res.insert({r.center, r.width, r.height, res.cache_entry->index});
  return res;
 }

 // HDF5
 friend std::string get_triqs_hdf5_data_scheme(configuration const&) { return "SomConfiguration"; }
 friend void h5_write(h5::group gr, std::string const& name, configuration const& c) {
  using triqs::arrays::array;
  array<double,2> data(c.size(),3);
  for(int i = 0; i < c.size(); ++i)
   data(i,range()) = array<double,1>{c[i].center, c[i].width, c[i].height};
  h5_write(gr, name, data);
  auto ds = gr.open_dataset(name);
  h5_write_attribute(ds, "TRIQS_HDF5_data_scheme", get_triqs_hdf5_data_scheme(c));
 }
 friend void h5_read(h5::group gr, std::string const& name, configuration& c, cache_index & ci) {
  using triqs::arrays::array;
  array<double,2> data;
  h5_read(gr, name, data);
  for(int i = 0; i < first_dim(data); ++i) c.insert({data(i,0), data(i,1), data(i,2), ci});
 }
};

}
