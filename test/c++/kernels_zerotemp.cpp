/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include <som/kernels/zerotemp_imfreq.hpp>
#include <som/kernels/zerotemp_imtime.hpp>
#include <som/kernels/zerotemp_legendre.hpp>

#include "./test_kernel.hpp"

class ZeroTemp_test : public ::testing::Test {
protected:
  cache_index ci;
};

TEST_F(ZeroTemp_test, imtime) {
  test_kernel<kernel<ZeroTemp, imtime>>("zerotemp_imtime.h5", ci, 1e-10);
}
TEST_F(ZeroTemp_test, imfreq) {
  test_kernel<kernel<ZeroTemp, imfreq>>("zerotemp_imfreq.h5", ci, 1e-10);
}
TEST_F(ZeroTemp_test, legendre) {
  test_kernel<kernel<ZeroTemp, legendre>>("zerotemp_legendre.h5", ci, 1e-10);
}
