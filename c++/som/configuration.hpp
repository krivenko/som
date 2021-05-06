/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2021 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <initializer_list>
#include <iostream>
#include <vector>

#include <mpi/mpi.hpp>
#include <mpi/vector.hpp>
#include <h5/h5.hpp>

#include "cache_index.hpp"
#include "rectangle.hpp"

namespace som {

struct configuration {

  // Rectangles in this configuration
  std::vector<rectangle> rects;

  static const int default_max_rects = 70;

  // Reference to the cache index and id within the cache.
  // Every new configuration object, including copies, acquire a new
  // cache entry descriptor, which is then released in the destructor.
  // Assignments and compound operations *=, += will invalidate the
  // cache entry. One should manually call kernel_base::cache_* methods
  // to update the invalidated entries without a full recomputation
  cache_index& ci;
  int cache_id;

private:
  // Configurations are supposed to be modified
  // only by config_update and solution_worker
  friend struct config_update;
  template <typename KernelType> friend class solution_worker;

  // Insert a new rectangle
  rectangle& insert(rectangle const& r);

  // Remove a rectangle by index
  void remove(int index);

  // Replace a rectangle
  void replace(int index, rectangle const& r);

  // Remove all rectangles, only for internal use in h5_read
  void clear();

public:
  explicit configuration(cache_index& ci,
                         int reserved_rects = default_max_rects);
  configuration(std::initializer_list<rectangle> const& l, cache_index& ci);
  configuration(configuration const& c);
  configuration(configuration&& c) noexcept;
  configuration& operator=(configuration const& c);
  configuration& operator=(configuration&& c) noexcept;
  ~configuration();

  // Number of rectangles
  [[nodiscard]] int size() const { return rects.size(); }
  // Maximum number of rectangles
  [[nodiscard]] int max_size() const { return rects.capacity(); }
  // Norm of this configuration
  [[nodiscard]] double norm() const;

  // Access a rectangle by index
  rectangle const& operator[](int index) const { return rects[index]; }

  // Equality
  bool operator==(configuration const& c) const { return rects == c.rects; }
  bool operator!=(configuration const& c) const { return !operator==(c); }

  // Sum of configurations: all rectangles from both of them
  configuration& operator+=(configuration const& c);
  configuration operator+(configuration const& c) const;

  // Multiply configuration by a positive scalar
  // Heights of all rectangles are scaled
  configuration& operator*=(double alpha);
  friend configuration operator*(configuration const& c, double alpha);
  friend configuration operator*(double alpha, configuration const& c);

  // Normalize configuration to have a total area of norm
  void normalize(double norm = 1.0);

  // constant iterator
  using const_iterator = std::vector<rectangle>::const_iterator;
  [[nodiscard]] const_iterator begin() const { return rects.begin(); }
  [[nodiscard]] const_iterator end() const { return rects.end(); }
  [[nodiscard]] const_iterator cbegin() const { return rects.cbegin(); }
  [[nodiscard]] const_iterator cend() const { return rects.cend(); }

  // stream insertion
  friend std::ostream& operator<<(std::ostream& os, configuration const& c);

  // MPI reduce
  friend configuration
  mpi_reduce(configuration const& c, mpi::communicator comm = {}, int root = 0,
             bool all = false,
             MPI_Op = MPI_SUM // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  ) {
    if(comm.size() == 1) return c;
    configuration res(c.ci);

    std::vector<rectangle::pod_t> pod_rects(std::begin(c), std::end(c));
    pod_rects = mpi_gather(pod_rects, comm, root, all);

    for(auto const& r : pod_rects)
      res.insert({r.center, r.width, r.height, res.ci});
    return res;
  }

  // HDF5
  friend std::string get_triqs_hdf5_data_scheme(configuration const&) {
    return "SomConfiguration";
  }
  friend void h5_write(h5::group gr, std::string const& name,
                       configuration const& c);
  friend void h5_read(h5::group gr, std::string const& name,
                      configuration& c, cache_index& ci_);
};

} // namespace som
