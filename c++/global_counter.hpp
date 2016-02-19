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

#include <triqs/mpi/base.hpp>

namespace som {

using triqs::mpi::communicator;

// Simple counter syncronized between all processes in a communicator
class global_counter {
 long * p_count;
 MPI_Win win;
 bool is_rank0;

public:

global_counter(communicator c, long count = 0) : is_rank0(c.rank() == 0) {
 MPI_Alloc_mem(sizeof(long),MPI_INFO_NULL,&p_count);
 *p_count = count;
 if(is_rank0)
  MPI_Win_create(p_count, sizeof(long), sizeof(long), MPI_INFO_NULL, c.get(), &win);
 else
  MPI_Win_create(p_count, 0, sizeof(long), MPI_INFO_NULL, c.get(), &win);
}

~global_counter() {
 MPI_Win_free(&win);
 MPI_Free_mem(p_count);
}

global_counter(global_counter const&) = delete;
global_counter& operator=(global_counter const&) = delete;
inline global_counter(global_counter && gc) {
 MPI_Free_mem(p_count);
 p_count = gc.p_count; win = gc.win; is_rank0 = gc.is_rank0;
}
inline global_counter& operator=(global_counter && gc) {
 MPI_Free_mem(p_count);
 p_count = gc.p_count; win = gc.win; is_rank0 = gc.is_rank0;
 return *this;
}

global_counter& operator=(long count) {
 MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
 MPI_Put(&count, 1, MPI_LONG, 0, 0, 1, MPI_LONG, win);
 MPI_Win_unlock(0, win);
 return *this;
}

inline global_counter& operator+=(long count) {
 MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
 if(!is_rank0) {
  MPI_Accumulate(&count, 1, MPI_LONG, 0, 0, 1, MPI_LONG, MPI_SUM, win);
 } else {
  (*p_count) += count;
 }
 MPI_Win_unlock(0, win);
 return *this;
}
inline global_counter& operator-=(long count) { return operator+=(-count); }
inline global_counter& operator++() { return operator+=(1); }
inline global_counter& operator--() { return operator+=(-1); }

inline explicit operator long() {
 if(!is_rank0) {
  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
  MPI_Get(p_count, 1, MPI_LONG, 0, 0, 1, MPI_LONG, win);
  MPI_Win_unlock(0, win);
 }
 return *p_count;
}

};

}
