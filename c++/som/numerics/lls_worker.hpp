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

#include <nda/nda.hpp>

namespace som {

// Solution of a linear least squares problem with a possibly rank-deficient
// matrix A and a right hand side b.
class lls_worker {

  int M;
  int N;
  nda::array<double, 1> b;
  nda::array<double, 1> s;
  nda::array<double, 1> work;

public:
  lls_worker(int n_eqs, int n_vars);

  int operator()(nda::matrix<double, nda::F_layout>& A,
                 nda::vector<double> const& b,
                 nda::vector<double>& x,
                 double rcond = -1);

  [[nodiscard]] inline nda::array<double, 1> const& singular_values() const {
    return s;
  }
};

} // namespace som
