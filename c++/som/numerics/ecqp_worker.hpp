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
#pragma once

#include <complex>
#include <type_traits>

#include <nda/nda.hpp>

namespace som {

/// This class wraps a few LAPACK calls to solve the Equality Constrained
/// Quadratic Programming (ECQP) problem:
///
/// Minimize
///     x^T \hat Q x - f^\dagger x - x^T f
/// w.r.t. real N-dimensional vector x subject to constraint
///     \hat L x = d.
///
/// Q is a real symmetric positive definite NxN matrix. f is a real vector of
/// size N, L is a real matrix of size PxN, and d is a real vector of size P.

class ecqp_worker {

  int N;
  int P;
  int gglse_ldb;

  nda::matrix<double, nda::F_layout> A;
  nda::vector<double> b;
  nda::matrix<double, nda::F_layout> L_;
  nda::vector<double> d_;
  nda::vector<double> work;

public:
  ecqp_worker(int N, int N_constraints);

  double operator()(nda::matrix_const_view<double, nda::F_layout> Q,
                    nda::vector_const_view<double> f,
                    nda::matrix_const_view<double, nda::F_layout> L,
                    nda::vector_const_view<double> d,
                    nda::vector_view<double> x);
};

} // namespace som
