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
#include <cmath>
#include <complex>
#include <numbers>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <som/numerics/dilog.hpp>

using namespace std::complex_literals;
using som::dilog;
using std::log;
using std::numbers::ln2;
using std::numbers::pi;

TEST(dilog_test, real) {
  EXPECT_CLOSE(-7.323953199000 - 14.467568824831i, dilog(100.0));
  EXPECT_CLOSE(0.5363012873579 - 7.2337844124155i, dilog(10.0));
  EXPECT_CLOSE(1.7837191612666 - 5.0561983221119i, dilog(5.0));
  EXPECT_CLOSE(2.4674010977725 - 2.1777431660094i, dilog(2.0001));
  EXPECT_CLOSE(pi * pi / 4 - 1i * pi * ln2, dilog(2.0));
  EXPECT_CLOSE(1.7423032794139 - 0.0622117884355i, dilog(1.02));
  EXPECT_CLOSE(1.7007321443240 - 0.0312598863091i, dilog(1.01));
  EXPECT_CLOSE(1.6459550523369 - 0.0003141435584i, dilog(1.0001));
  EXPECT_CLOSE(pi * pi / 6, dilog(1));
  EXPECT_CLOSE(0.9784693929303, dilog(0.75));
  EXPECT_CLOSE(pi * pi / 12 - 0.5 * ln2 * ln2, dilog(0.5));
  EXPECT_CLOSE(0.3261295100755, dilog(0.3));
  EXPECT_CLOSE(0.2676526390827, dilog(0.25));
  EXPECT_CLOSE(0.1026177910994, dilog(0.1));
  EXPECT_CLOSE(0.0001000025001, dilog(0.0001));
  EXPECT_CLOSE(0, dilog(0));
  EXPECT_CLOSE(-0.000099997500260, dilog(-0.0001));
  EXPECT_CLOSE(-0.097605235229322, dilog(-0.1));
  EXPECT_CLOSE(-pi * pi / 12, dilog(-1));
  EXPECT_CLOSE(-4.198277886858104, dilog(-10));
  EXPECT_CLOSE(-12.23875517731494, dilog(-100));
}
