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
#include <som/kernels/bosoncorr_imfreq.hpp>
#include <som/kernels/bosoncorr_imtime.hpp>
#include <som/kernels/bosoncorr_legendre.hpp>

#include "./test_kernel.hpp"

class BosonCorr_test : public ::testing::Test {
protected:
  cache_index ci;
};

TEST_F(BosonCorr_test, imtime) {
  test_kernel<kernel<BosonCorr, imtime>>("bosoncorr_imtime.h5", ci, 1e-10);
}
TEST_F(BosonCorr_test, imfreq) {
  test_kernel<kernel<BosonCorr, imfreq>>("bosoncorr_imfreq.h5", ci, 1e-10);
}
TEST_F(BosonCorr_test, legendre) {
  test_kernel<kernel<BosonCorr, legendre>>("bosoncorr_legendre.h5", ci, 1e-10);
}
