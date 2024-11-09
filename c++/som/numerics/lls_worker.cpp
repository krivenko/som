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

#include <cassert>
#include <cmath>

#include <nda/lapack/interface/lapack.h>

#include <triqs/utility/exceptions.hpp>

#include "lls_worker.hpp"

namespace som {
namespace lapack::f77 {

  void gelss(int M,
             int N,
             int NRHS,
             double *A,
             int LDA,
             double *B,
             int LDB,
             double *S,
             double RCOND,
             int &RANK,
             double *WORK,
             int LWORK,
             int &INFO) {
    LAPACK_dgelss(
        &M, &N, &NRHS, A, &LDA, B, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);
  }

} // namespace lapack::f77

lls_worker::lls_worker(int n_eqs, int n_vars) : M(n_eqs), N(n_vars) {
  if(n_eqs <= 0)
    TRIQS_RUNTIME_ERROR
        << "lls_worker: The number of equations must be positive, got "
        << n_eqs;
  if(n_vars <= 0)
    TRIQS_RUNTIME_ERROR
        << "lls_worker: The number of variables must be positive, got "
        << n_vars;

  b.resize(std::max(n_eqs, n_vars));
  s.resize(std::min(n_eqs, n_vars));

  // Query optimal value of LWORK
  int rank = 0;
  double opt_lwork = 0;
  int info = 0;
  lapack::f77::gelss(M,
                     N,
                     1,
                     nullptr,
                     n_eqs,
                     nullptr,
                     int(b.size()),
                     nullptr,
                     -1,
                     rank,
                     &opt_lwork,
                     -1,
                     info);
  if(info == 0)
    work.resize(long(opt_lwork));
  else
    TRIQS_RUNTIME_ERROR
        << "lls_worker: Could not calculate the optimal size of the WORK array";
}

int lls_worker::operator()(nda::matrix<double, nda::F_layout> &A,
                           nda::vector<double> const &rhs,
                           nda::vector<double> &x,
                           double rcond) {
  assert(first_dim(A) == M);
  assert(second_dim(A) == N);
  assert(rhs.size() == M);
  assert(x.size() == N);

  b(nda::range(M)) = rhs;
  int rank = 0;
  int info = 0;
  lapack::f77::gelss(M,
                     N,
                     1,
                     A.data(),
                     M,
                     b.data(),
                     int(b.size()),
                     s.data(),
                     rcond,
                     rank,
                     work.data(),
                     int(work.size()),
                     info);
  if(info != 0) TRIQS_RUNTIME_ERROR << "lls_worker: SVD failed to converge";
  x() = b(nda::range(N));
  return rank;
}

} // namespace som
