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

// MCS (Mellor-Crummey, Scott) mutex
// Used to establish critical sections for concurrent MPI RMA operations
//
// Implementation based on
// T. Hoefler et al., 2013. Remote Memory Access Programming in MPI-3.
// ACM Trans. Parallel Comput. 1, 1, Article 1 (March 2013), 29 pages.
// Requires MPI-3.0
class mpi_mutex {

 int * data;        // Pointer to (locked, next) pair
 MPI_Win win;       // RMA window associated with data
 communicator comm; // Communicator
 int tail_rank;     // Rank of the process in comm that holds the tail pointer
 long id;           // ID of this mutex

 static constexpr int mutex_released_tag_base = 1024;
 static long get_id(bool next) {
  static long id_ = 0;
  return next ? id_++ : id_--;
 }

 static constexpr size_t locked_disp = 0; // Displacement of the 'locked' flag in the data buffer
 static constexpr size_t next_disp = 1;   // Displacement of the 'next' field in the data buffer

public:

 mpi_mutex(communicator const& c, int tail_rank) : comm(c), tail_rank(tail_rank), id(get_id(true)) {
  // Allocate memory for data
  MPI_Win_allocate(2 * sizeof(int), sizeof(int), MPI_INFO_NULL, comm.get(), &data, &win);
  // Enable shared access to all processes in win
  MPI_Win_lock_all(0, win);

  // Initialize data
  data[0] = MPI_PROC_NULL;
  data[1] = MPI_PROC_NULL;
  MPI_Win_sync(win);

  comm.barrier();
 }

 ~mpi_mutex() {
  comm.barrier();
  get_id(false);

  // End RMA access epoch on win
  MPI_Win_unlock_all(win);
  // Free memory
  MPI_Win_free(&win);
 }

 mpi_mutex(mpi_mutex const&) = delete;
 mpi_mutex(mpi_mutex &&) = delete;
 mpi_mutex & operator=(mpi_mutex const&) = delete;
 mpi_mutex & operator=(mpi_mutex &&) = delete;

 void lock() {
  int rank = comm.rank();

  // This store is safe, since it cannot happen concurrently with a remote write
  data[locked_disp] = MPI_PROC_NULL;
  MPI_Win_sync(win);

  int prev;
  MPI_Fetch_and_op(&rank, &prev, MPI_INT, tail_rank, next_disp, MPI_REPLACE, win);
  MPI_Win_flush(tail_rank, win);

  // If there was a previous tail, update their next pointer and wait for
  // notification.  Otherwise, the mutex was successfully acquired.
  if (prev != MPI_PROC_NULL) {
   // Wait for notification

   MPI_Accumulate(&rank, 1, MPI_INT, prev, locked_disp, 1, MPI_INT, MPI_REPLACE, win);
   MPI_Win_flush(prev, win);

   MPI_Recv(nullptr, 0, MPI_BYTE, prev, mutex_released_tag_base + id, comm.get(), MPI_STATUS_IGNORE);
  }
 }

 void unlock() {
  MPI_Win_sync(win);

  int rank = comm.rank(), next;

  // Read my next pointer.
  MPI_Fetch_and_op(NULL, &next, MPI_INT, comm.rank(), locked_disp, MPI_NO_OP, win);
  MPI_Win_flush(rank, win);

  if(next == MPI_PROC_NULL) {
   int tail;
   int nil = MPI_PROC_NULL;

   // Check if we are the at the tail of the lock queue. If so, we are done, if not, we need to send notification.
   MPI_Compare_and_swap(&nil, &rank, &tail, MPI_INT, tail_rank, next_disp, win);
   MPI_Win_flush(tail_rank, win);

   if(tail != rank) {
    for (;;) {
     MPI_Fetch_and_op(nullptr, &next, MPI_INT, comm.rank(), locked_disp, MPI_NO_OP, win);
     MPI_Win_flush(rank, win);
     if(next != MPI_PROC_NULL) break;
    }
   }
  }

  // Notify the next waiting process
  if (next != MPI_PROC_NULL)
   MPI_Send(nullptr, 0, MPI_BYTE, next, mutex_released_tag_base + id, comm.get());
 }
};

}
