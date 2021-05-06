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

#include <array>
#include <cmath>

#include <triqs/utility/numeric_ops.hpp>

#include "rectangle.hpp"

namespace som {

using triqs::arrays::vector;

bool rectangle::operator==(rectangle const& r) const {
  using triqs::utility::is_zero;
  return is_zero(center - r.center, center_equal_tol) &&
         is_zero(width - r.width, width_equal_tol) &&
         is_zero(height - r.height, height_equal_tol);
}

double rectangle::operator()(double x) const {
  return (x >= center - width / 2 && x <= center + width / 2) ? height : 0;
}

std::complex<double> rectangle::hilbert_transform(std::complex<double> z,
                                                  bool multiply_by_e) const {
  if(multiply_by_e)
    // -h \int_{c-w/2}^{c+w/2} d\epsilon'
    //   \frac{\epsilon'}{\epsilon' - \epsilon - i0}
    return -height * (width + z * std::log((center + width / 2 - z) /
                                           (center - width / 2 - z)));
  else
    // -h \int_{c-w/2}^{c+w/2} d\epsilon' \frac{1}{\epsilon' - \epsilon - i0}
    return -height *
           std::log((center + width / 2 - z) / (center - width / 2 - z));
}

vector<double> rectangle::tail_coefficients(long order_min, long order_max,
                                            bool multiply_by_e) const {
  vector<double> data(order_max - order_min + 1);
  double e1 = center - width / 2, e2 = center + width / 2;
  double e1n, e2n;
  int denom_shift;
  if(multiply_by_e) {
    e1n = e1;
    e2n = e2;
    denom_shift = 1;
  } else {
    e1n = 1.0;
    e2n = 1.0;
    denom_shift = 0;
  }
  for(long n = order_min; n <= order_max; ++n) {
    if(n < 1) {
      data(n - order_min) = 0;
      continue;
    }
    e1n *= e1;
    e2n *= e2;
    data(n - order_min) = height * (e2n - e1n) / double(n + denom_shift);
  }
  return data;
}

// Output stream insertion
std::ostream& operator<<(std::ostream& os, rectangle const& r) {
  os << "(c:" << r.center << ", w:" << r.width << ", h:" << r.height << ")";
  return os;
}

} // namespace som

namespace mpi {

MPI_Datatype mpi_type<som::rectangle::pod_t>::get() noexcept {
  static bool type_committed = false;
  static MPI_Datatype dt;
  if(!type_committed) {
    std::array<int, 3> blocklengths = {1, 1, 1};
    std::array<MPI_Aint, 3> displacements = {0, sizeof(double),
                                             2 * sizeof(double)};
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    std::array<MPI_Datatype, 3> types = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(3, blocklengths.data(), displacements.data(),
                           types.data(), &dt);
    MPI_Type_commit(&dt);
    type_committed = true;
  }
  return dt;
}

} // namespace mpi
