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
#include <som/kernels/fermiongfsymm_imfreq.hpp>
#include <som/kernels/fermiongfsymm_imtime.hpp>
#include <som/kernels/fermiongfsymm_legendre.hpp>

#include "./test_kernel.hpp"

class FermionGfSymm_test : public ::testing::Test {
protected:
  cache_index ci;
};

TEST_F(FermionGfSymm_test, imtime) {
  test_kernel<kernel<FermionGfSymm, imtime>>(
      "fermiongfsymm_imtime.h5", ci, 1e-10);
}
TEST_F(FermionGfSymm_test, imfreq) {
  test_kernel<kernel<FermionGfSymm, imfreq>>(
      "fermiongfsymm_imfreq.h5", ci, 1e-10);
}
TEST_F(FermionGfSymm_test, legendre) {
  test_kernel<kernel<FermionGfSymm, legendre>>(
      "fermiongfsymm_legendre.h5", ci, 1e-10);
}
