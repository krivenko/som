/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko
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

#include <numeric>
#include <utility>

#include <mpi/vector.hpp>

#include "global_index_map.hpp"

namespace som {

std::vector<std::size_t> init_range_size(mpi::communicator const& comm,
                                         std::size_t local_block_size) {
  std::vector<std::size_t> tmp{local_block_size};
  return mpi::all_gather(std::move(tmp), comm);
}

global_index_map::global_index_map(mpi::communicator const& comm,
                                   std::size_t local_block_size)
   : range_start_(comm.size(), 0)
   , range_size_(init_range_size(comm, local_block_size))
   , size_(std::accumulate(range_size_.begin(), range_size_.end(), 0)) {
  std::partial_sum(
      range_size_.begin(), range_size_.end() - 1, range_start_.begin() + 1);
}

} // namespace som
