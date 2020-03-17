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
#include <cmath>
#include <triqs/test_tools/arrays.hpp>

#include <som/numerics/dilog.hpp>

using som::dilog;
using std::log;

TEST(dilog, real) {
  EXPECT_CLOSE(-7.323953199000 - 14.467568824831_j, dilog(100.0));
  EXPECT_CLOSE(0.5363012873579 - 7.2337844124155_j, dilog(10.0));
  EXPECT_CLOSE(1.7837191612666 - 5.0561983221119_j, dilog(5.0));
  EXPECT_CLOSE(2.4674010977725 - 2.1777431660094_j, dilog(2.0001));
  EXPECT_CLOSE(M_PI * M_PI / 4 - 1_j * M_PI * log(2), dilog(2.0));
  EXPECT_CLOSE(1.7423032794139 - 0.0622117884355_j, dilog(1.02));
  EXPECT_CLOSE(1.7007321443240 - 0.0312598863091_j, dilog(1.01));
  EXPECT_CLOSE(1.6459550523369 - 0.0003141435584_j, dilog(1.0001));
  EXPECT_CLOSE(M_PI * M_PI / 6, dilog(1));
  EXPECT_CLOSE(0.9784693929303, dilog(0.75));
  EXPECT_CLOSE(M_PI * M_PI / 12 - 0.5 * log(2) * log(2), dilog(0.5));
  EXPECT_CLOSE(0.3261295100755, dilog(0.3));
  EXPECT_CLOSE(0.2676526390827, dilog(0.25));
  EXPECT_CLOSE(0.1026177910994, dilog(0.1));
  EXPECT_CLOSE(0.0001000025001, dilog(0.0001));
  EXPECT_CLOSE(0, dilog(0));
  EXPECT_CLOSE(-0.000099997500260, dilog(-0.0001));
  EXPECT_CLOSE(-0.097605235229322, dilog(-0.1));
  EXPECT_CLOSE(-M_PI * M_PI / 12, dilog(-1));
  EXPECT_CLOSE(-4.198277886858104, dilog(-10));
  EXPECT_CLOSE(-12.23875517731494, dilog(-100));
}

MAKE_MAIN;
