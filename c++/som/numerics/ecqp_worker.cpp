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

#include <cassert>
#include <cmath>

#include <nda/blas/interface/cblas_f77.h>
#include <nda/lapack/interface/lapack.h>

#include <triqs/utility/exceptions.hpp>

#include "ecqp_worker.hpp"

namespace som {
namespace blas::f77 {

  //
  // Solution of a linear system with a triangular matrix
  //

  void trsm(char SIZE,
            char UPLO,
            char TRANSA,
            char DIAG,
            int M,
            int N,
            double ALPHA,
            const double *A,
            int LDA,
            double *B,
            int LDB) {
    F77_dtrsm(&SIZE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, A, &LDA, B, &LDB);
  }

} // namespace blas::f77

namespace lapack::f77 {

  //
  // Cholesky factorization of a positive definite matrix
  //

  void potrf(char UPLO, int N, double *A, int LDA, int &INFO) {
    LAPACK_dpotrf(&UPLO, &N, A, &LDA, &INFO);
  }

  //
  // Constrained least squares
  //

  void gglse(int M,
             int N,
             int P,
             double *A,
             int LDA,
             double *B,
             int LDB,
             double *C,
             double *D,
             double *X,
             double *WORK,
             int LWORK,
             int &INFO) {
    LAPACK_dgglse(&M, &N, &P, A, &LDA, B, &LDB, C, D, X, WORK, &LWORK, &INFO);
  }

} // namespace lapack::f77

ecqp_worker::ecqp_worker(int N, int N_constraints)
   : N(N), P(N_constraints), gglse_ldb(std::max(1, P)) {
  if(N <= 0)
    TRIQS_RUNTIME_ERROR
        << "glse_worker: The number of variables must be positive, got " << N;
  if(P < 0)
    TRIQS_RUNTIME_ERROR
        << "ecqp_worker: The number of constraints must be non-negative, got "
        << P;
  if(P > N)
    TRIQS_RUNTIME_ERROR << "ecqp_worker: The number of constraints cannot "
                           "exceed the number of variables ("
                        << N << " vs " << P << ")";
  A.resize(N, N);
  b.resize(N);
  L_.resize(gglse_ldb, N);
  d_.resize(P);

  // Query optimal value of LWORK
  double opt_lwork = 0;
  int info = 0;
  lapack::f77::gglse(N,
                     N,
                     P,
                     nullptr,
                     N,
                     nullptr,
                     gglse_ldb,
                     nullptr,
                     nullptr,
                     nullptr,
                     &opt_lwork,
                     -1,
                     info);
  if(info == 0)
    work.resize(long(opt_lwork));
  else
    TRIQS_RUNTIME_ERROR << "ecqp_worker: Could not calculate the optimal size "
                           "of the WORK array";
}

double ecqp_worker::operator()(nda::matrix_const_view<double, nda::F_layout> Q,
                               nda::vector_const_view<double> f,
                               nda::matrix_const_view<double, nda::F_layout> L,
                               nda::vector_const_view<double> d,
                               nda::vector_view<double> x) {
  assert(first_dim(Q) == N);
  assert(second_dim(Q) == N);
  assert(f.size() == N);
  assert(first_dim(L) == P);
  assert(second_dim(L) == N);
  assert(d.size() == P);
  assert(x.size() == N);

  // Perform Cholesky factorization of Q
  A() = Q;
  int info = 0;
  lapack::f77::potrf('U', N, A.data(), N, info);
  if(info != 0) {
    if(info > 0)
      TRIQS_RUNTIME_ERROR << "ecqp_worker: The leading minor of order " << info
                          << " is not positive definite";
    else
      TRIQS_RUNTIME_ERROR << "ecqp_worker: The " << (-info)
                          << "-th argument of ?potrf had an illegal value";
  }

  // Fill the strictly lower triangular part of A with zeros
  for(auto i : nda::range(1, N))
    for(auto j : nda::range(0, i)) A(i, j) = 0;

  // Compute RHS of the least squares problem
  b() = f;
  blas::f77::trsm('L', 'U', 'T', 'N', N, 1, 1.0, A.data(), N, b.data(), N);
  auto b2 = sum(abs2(b));

  // Solve the constrained least squares
  if(P > 0) {
    L_() = L;
    d_() = d;
  }
  info = 0;
  lapack::f77::gglse(N,
                     N,
                     P,
                     A.data(),
                     N,
                     L_.data(),
                     gglse_ldb,
                     b.data(),
                     d_.data(),
                     x.data(),
                     work.data(),
                     int(work.size()),
                     info);
  if(info != 0) {
    switch(info) {
      case 1:
        TRIQS_RUNTIME_ERROR
            << "ecqp_worker: The constraints matrix is rank-deficient";
        break;
      case 2:
        TRIQS_RUNTIME_ERROR
            << "ecqp_worker: The combined LS/constraints matrix "
               "does not have the full column rank";
        break;
      default:
        TRIQS_RUNTIME_ERROR << "ecqp_worker: The " << (-info)
                            << "-th argument of ?gglse() had an illegal value";
    }
  }

  // Return the value of the objective function
  // x^T \hat Q x - f^\dagger x - x^T f at the minimum.
  // It equals | A*x - b |^2 - b^2.
  return sum(abs2(b(nda::range(N - P, N)))) - b2;
}

} // namespace som
