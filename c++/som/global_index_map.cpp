/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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

global_index_map::global_index_map(mpi::communicator const& comm,
                                   int local_block_size)
   : range_start_(comm.size(), 0), range_size_(comm.size()), size_(0) {
  std::vector<int> tmp{local_block_size};
  range_size_ = mpi::all_gather(std::move(tmp), comm);
  std::partial_sum(range_size_.begin(),
                   range_size_.end() - 1,
                   range_start_.begin() + 1);
  size_ = std::accumulate(range_size_.begin(), range_size_.end(), 0);
}

} // namespace som
