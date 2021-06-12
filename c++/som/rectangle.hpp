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

#include <complex>
#include <ostream>
#include <utility>

#include <mpi/mpi.hpp>
#include <triqs/arrays/vector.hpp>

#include "cache_index.hpp"

namespace som {

struct rectangle {

  double center; // center
  double width;  // width
  double height; // height

  // Pointer to the associated cache entry (can be null)
  //
  // Rectangles are immutable objects, therefore they acquire a new cache entry
  // descriptor only upon construction. Copies inherit the same ID, and
  // increase the refcount to it. The refcount is decreased by one on
  // destruction. Cache entry bound to an existing rectangle is never
  // invalidated, once the LHS is computed.
  cache_entry_ptr cache_ptr;

public:

  // Construct a new rectangle unbound from a LHS cache
  inline rectangle(double center, double width, double height)
     : center(center)
     , width(width)
     , height(height)
     {}

  // Construct a new rectangle and bind it to a new cache entry
  inline rectangle(double center, double width, double height, cache_index& ci)
     : center(center)
     , width(width)
     , height(height)
     , cache_ptr(ci)
     {}

  rectangle(rectangle const& r) = default;
  rectangle(rectangle&& r) noexcept = default;

  inline rectangle& operator=(rectangle const& r) {
    center = r.center;
    width = r.width;
    height = r.height;
    cache_ptr = r.cache_ptr;
    return *this;
  }
  inline rectangle& operator=(rectangle&& r) noexcept {
    using std::swap;
    swap(center, r.center);
    swap(width, r.width);
    swap(height, r.height);
    cache_ptr = std::move(r.cache_ptr);
    return *this;
  }

  static constexpr double center_equal_tol = 1e-8;
  static constexpr double width_equal_tol = 1e-8;
  static constexpr double height_equal_tol = 1e-8;

  bool operator==(rectangle const& r) const;
  bool operator!=(rectangle const& r) const { return !operator==(r); }

  [[nodiscard]] double norm() const { return width * height; }
  double operator()(double x) const;

  [[nodiscard]] std::complex<double>
  hilbert_transform(std::complex<double> z, bool multiply_by_e = false) const;
  [[nodiscard]] triqs::arrays::vector<double>
  tail_coefficients(long order_min, long order_max,
                    bool multiply_by_e = false) const;

  // Output stream insertion
  friend std::ostream& operator<<(std::ostream& os, rectangle const& r);

  // Multiplication by scalar
  friend rectangle operator*(rectangle const& r, double alpha) {
    return {r.center, r.width, r.height * alpha, r.cache_ptr.get_ci()};
  }
  friend rectangle operator*(double alpha, rectangle const& r) {
    return r * alpha;
  }

  // POD version of the rectangle, used in MPI operations
  struct pod_t {
    double center;
    double width;
    double height;
    pod_t() = default;
    explicit pod_t(rectangle const& r)
       : center(r.center), width(r.width), height(r.height) {}
  };

  // Convert to tuple (center,width,height)
  operator std::tuple<double, double, double>() const {
    return std::make_tuple(center, width, height);
  }
};

} // namespace som

namespace mpi {

template <> struct mpi_type<som::rectangle::pod_t> {
  static MPI_Datatype get() noexcept;
};

} // namespace mpi
