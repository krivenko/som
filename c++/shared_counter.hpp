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
#include "mpi_mutex.hpp"

namespace som {

using triqs::mpi::communicator;

// Simple counter syncronized between all processes in a communicator
//
// Many MPI implementations support passive target synchronization rather poorly.
// That is why calls to MPI_Win_flush/MPI_Win_lock/MPI_Win_unlock can block
// until an MPI function is called on the target side.
class shared_counter {
 long * p_count;
 MPI_Win win;
 communicator comm;
 mpi_mutex mutex;

public:

shared_counter(communicator const& c, long count = 0) : comm(c), mutex(comm, 0) {
 MPI_Alloc_mem(sizeof(long), MPI_INFO_NULL, &p_count);
 if(c.rank() == 0)
  MPI_Win_create(p_count, sizeof(long), sizeof(long), MPI_INFO_NULL, c.get(), &win);
 else
  MPI_Win_create(p_count, 0, sizeof(long), MPI_INFO_NULL, c.get(), &win);

 // Lock the window for shared access, we will use a mutex to implement critical sections
 MPI_Win_lock(MPI_LOCK_SHARED, 0, MPI_MODE_NOCHECK, win);
 if(c.rank() == 0) {
  MPI_Put(&count, 1, MPI_LONG, 0, 0, 1, MPI_LONG, win);
  MPI_Win_flush(0, win);
 }
 comm.barrier();
}

~shared_counter() {
 comm.barrier();
 MPI_Win_unlock(0, win);
 MPI_Win_free(&win);
 MPI_Free_mem(p_count);
}

shared_counter(shared_counter const&) = delete;
shared_counter& operator=(shared_counter const&) = delete;
inline shared_counter(shared_counter && gc) = delete;
inline shared_counter& operator=(shared_counter && gc) = delete;

shared_counter& operator=(long count) {
 mutex.lock();
 MPI_Win_sync(win);
 MPI_Put(&count, 1, MPI_LONG, 0, 0, 1, MPI_LONG, win);
 MPI_Win_flush(0, win);
 mutex.unlock();
 return *this;
}

inline shared_counter& operator+=(long count) {
 mutex.lock();
 MPI_Win_sync(win);
 MPI_Accumulate(&count, 1, MPI_LONG, 0, 0, 1, MPI_LONG, MPI_SUM, win);
 MPI_Win_flush(0, win);
 mutex.unlock();
 return *this;
}
inline shared_counter& operator-=(long count) { return operator+=(-count); }
inline shared_counter& operator++() { return operator+=(1); }
inline shared_counter& operator--() { return operator+=(-1); }

// Postfix changes of the counter
inline long postfix_add(long count) {
 mutex.lock();
 MPI_Win_sync(win);
 long old_val;
 MPI_Fetch_and_op(&count, &old_val, MPI_LONG, 0, 0, MPI_SUM, win);
 MPI_Win_flush(0, win);
 mutex.unlock();
 return old_val;
}
inline long operator++(int) { return postfix_add(1); }
inline long operator--(int) { return postfix_add(-1); }

inline explicit operator long() {
 mutex.lock();
 MPI_Win_sync(win);
 long val;
 MPI_Get(&val, 1, MPI_LONG, 0, 0, 1, MPI_LONG, win);
 MPI_Win_flush(0, win);
 mutex.unlock();
 return val;
}

};

}
