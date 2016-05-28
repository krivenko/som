/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016 by I. Krivenko
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
#include <triqs/test_tools/arrays.hpp>

#include "polynomial.hpp"

using som::polynomial;

TEST(polynomial, empty) {
 polynomial p;

 EXPECT_EQ(0,p(0));
 EXPECT_EQ(0,p(10));
 EXPECT_EQ(0,p(76.2));
}

TEST(polynomial, degree5) {
 polynomial p({1., 3., -2., 5., -11., 0.4});

 EXPECT_NEAR(-3.6,           p(1.0),  1e-10);
 EXPECT_NEAR(-65169,         p(10.0), 1e-10);
 EXPECT_NEAR(-144047.192472, p(12.7), 1e-10);
}

MAKE_MAIN;
