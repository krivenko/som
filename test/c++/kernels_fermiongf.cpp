/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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
#include "kernels/fermiongf_imtime.hpp"
#include "kernels/fermiongf_imfreq.hpp"
#include "kernels/fermiongf_legendre.hpp"

#include "./test_kernel.hpp"

cache_index ci;

TEST(FermionGf, imtime) {
 test_kernel<kernel<FermionGf,imtime>>("fermiongf_imtime.h5", ci, 1e-10);
}
TEST(FermionGf, imfreq) {
 test_kernel<kernel<FermionGf,imfreq>>("fermiongf_imfreq.h5", ci, 1e-10);
}
TEST(FermionGf, legendre) {
 test_kernel<kernel<FermionGf,legendre>>("fermiongf_legendre.h5", ci, 1e-10);
}

MAKE_MAIN;
