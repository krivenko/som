/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2025 Igor Krivenko
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

#include <iostream>

#include <gtest/gtest.h>

#include <som/global_index_map.hpp>

using namespace som;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  if(mpi::has_env) {
    mpi::environment env(argc, argv);
    std::cout << "MPI environment detected\n";
    return RUN_ALL_TESTS();
  } else
    return RUN_ALL_TESTS();
}

TEST(global_index_map, global_index_map) {
  mpi::communicator comm;
  int size = comm.size();

  global_index_map map(comm, comm.rank() + 1);

  ASSERT_EQ(map.size(), size * (1 + size) / 2);

  int global_index = 0;
  for(int rank = 0; rank < size; ++rank) {
    EXPECT_EQ(map.range_start(rank), rank * (1 + rank) / 2);
    EXPECT_EQ(map.range_size(rank), rank + 1);

    for(int local_index = 0; local_index < rank + 1; ++local_index) {
      EXPECT_EQ(map(rank, local_index), global_index);
      ++global_index;
    }
  }
}
