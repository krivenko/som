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
#include <triqs/test_tools/arrays.hpp>
#include <sstream>

#include "numerics/polynomial.hpp"

using som::polynomial;

TEST(polynomial, empty) {
 polynomial<> p;

 EXPECT_EQ(0,p(0));
 EXPECT_EQ(0,p(10));
 EXPECT_EQ(0,p(76.2));
}

TEST(polynomial, degree5) {
 polynomial<> p{1., 3., -2., 5., -11., 0.4};

 EXPECT_EQ(6, p.size());
 EXPECT_EQ(5, p.degree());
 EXPECT_NEAR(-3.6,           p(1.0),  1e-10);
 EXPECT_NEAR(-65169,         p(10.0), 1e-10);
 EXPECT_NEAR(-144047.192472, p(12.7), 1e-10);
}

TEST(polynomial, try_lower_degree) {
 polynomial<> p{1., 3., -2., 5., 0., 2.0, 0., 0.};

 EXPECT_EQ(7, p.degree());
 p.try_lower_degree();
 EXPECT_EQ(5, p.degree());

 p = polynomial<>{1e-20, 0., 0., 0.};
 EXPECT_EQ(4, p.size());
 p.try_lower_degree();
 EXPECT_EQ(0, p.size());
}

TEST(polynomial, print) {
 auto check_print = [](polynomial<> const& p, std::string const& ref) {
  std::stringstream ss;
  ss << p;
  EXPECT_EQ(ref, ss.str());
 };

 polynomial<> p_empty{};
 polynomial<> p0{2};
 polynomial<> p1{2, 3};
 polynomial<> p4{2, 7, 9, 11, 4};

 check_print(p_empty, "0");
 check_print(p0, "2");
 check_print(p1, "2 + 3*x");
 check_print(p4, "2 + 7*x + 9*x^2 + 11*x^3 + 4*x^4");
}

TEST(polynomial, arithmetics) {
 polynomial<> p_empty{};
 polynomial<> p1{1, 2, 3, 4};
 polynomial<> p2{1, 2, 3, 4, 5, 6};
 polynomial<> p3{7, 8, 9, 10};

 // Comparison
 EXPECT_EQ(p1, p1);
 EXPECT_NE(p1, p2);

 // Addition
 EXPECT_EQ(polynomial<>({2, 4, 6, 8, 5, 6}), p1+p2);
 EXPECT_EQ(polynomial<>({2, 4, 6, 8, 5, 6}), p2+p1);
 EXPECT_EQ(polynomial<>({8, 10, 12, 14}), p1+p3);

 // Subtraction
 EXPECT_EQ(polynomial<>({0, 0, 0, 0, -5, -6}), p1-p2);
 EXPECT_EQ(polynomial<>({0, 0, 0, 0, 5, 6}), p2-p1);
 EXPECT_EQ(polynomial<>({6, 6, 6, 6}), p3-p1);

 // Multiplication
 EXPECT_EQ(polynomial<>({}), p1*p_empty);
 EXPECT_EQ(polynomial<>({}), p_empty*p2);
 EXPECT_EQ(polynomial<>({1, 4, 10, 20, 30, 40, 43, 38, 24}), p1*p2);
 EXPECT_EQ(polynomial<>({}), 2.0*p_empty);
 EXPECT_EQ(polynomial<>({}), p_empty*2.0);
 EXPECT_EQ(polynomial<>({}), 2.0*p_empty);
 EXPECT_EQ(polynomial<>({2, 4, 6, 8}), 2.0*p1);
 EXPECT_EQ(polynomial<>({2, 4, 6, 8}), p1*2.0);
}

TEST(polynomial, derivative) {
 polynomial<> p_empty{};
 polynomial<> p0{4};
 polynomial<> p5{1, 2, 3, 4, 5, 6};

 EXPECT_EQ(polynomial<>({}), derivative(p_empty));
 EXPECT_EQ(polynomial<>({}), derivative(p0));
 EXPECT_EQ(polynomial<>({2, 6, 12, 20, 30}), derivative(p5));
}

TEST(polynomial, antiderivative) {
 polynomial<> p_empty{};
 polynomial<> p0{4};
 polynomial<> p5{1, 2, 3, 4, 5, 6};

 EXPECT_EQ(polynomial<>({}), antiderivative(p_empty));
 EXPECT_EQ(polynomial<>({0, 4}), antiderivative(p0));
 EXPECT_EQ(polynomial<>({0, 1, 1, 1, 1, 1, 1}), antiderivative(p5));
}

MAKE_MAIN;
