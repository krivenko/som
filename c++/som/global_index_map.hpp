/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <cassert>
#include <utility>
#include <vector>

#include <mpi/mpi.hpp>

namespace som {

// Mapping (rank, local index) -> global index
class global_index_map {

  // Parameters of global index ranges corresponding to MPI ranks.
  std::vector<std::size_t> range_start_;
  std::vector<std::size_t> range_size_;
  // Total number of global index values.
  std::size_t size_;

public:
  global_index_map(mpi::communicator const& comm, std::size_t local_block_size);

  [[nodiscard]] std::size_t range_start(int rank) const {
    return range_start_[rank];
  }
  [[nodiscard]] std::size_t range_size(int rank) const {
    return range_size_[rank];
  }
  [[nodiscard]] std::size_t size() const { return size_; }

  inline std::size_t operator()(int rank, std::size_t local_index) const {
    assert(local_index < range_size_[rank]);
    return range_start_[rank] + local_index;
  }
};

} // namespace som
