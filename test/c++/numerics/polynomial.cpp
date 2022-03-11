/*******************************************************************************
 *
 * SOM: Stochastic Optimization Method for Analytic Continuation
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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

#include <sstream>
#include <type_traits>

#include <nda/gtest_tools.hpp>

#include <som/numerics/polynomial.hpp>

using som::polynomial;

#define EXPECT_PRINT(X, Y)                                                     \
  {                                                                            \
    std::stringstream ss;                                                      \
    ss << Y;                                                                   \
    EXPECT_EQ(X, ss.str());                                                    \
  }

template <typename Polynomial> class polynomials_case : public ::testing::Test {

public:
  polynomials_case() = default;

  void test_empty() {
    Polynomial p;

    EXPECT_EQ(.0, p(0));
    EXPECT_EQ(.0, p(10));
    EXPECT_EQ(.0, p(76.2));
  }

  void test_degree5() {
    Polynomial p{1., 3., -2., 5., -11., 0.4};

    EXPECT_EQ(6, p.size());
    EXPECT_EQ(5, p.degree());

    if constexpr(std::is_same_v<typename Polynomial::coeff_type, double>) {
      EXPECT_NEAR(-3.6, p(1.0), 1e-10);
      EXPECT_NEAR(-65169, p(10.0), 1e-10);
      EXPECT_NEAR(-144047.192472, p(12.7), 1e-10);

    } else {
      EXPECT_COMPLEX_NEAR(-3.6, p(1.0), 1e-10);
      EXPECT_COMPLEX_NEAR(-65169, p(10.0), 1e-10);
      EXPECT_COMPLEX_NEAR(-144047.192472, p(12.7), 1e-10);
    }
  }

  void test_try_lower_degree() {
    Polynomial p{1., 3., -2., 5., 0., 2.0, 0., 0.};
/*
    if constexpr(Polynomial::auto_lower_degree)
      EXPECT_EQ(5, p.degree());
    else
      EXPECT_EQ(7, p.degree());
    p.try_lower_degree();
    EXPECT_EQ(5, p.degree());

    p = Polynomial{1e-20, 0., 0., 0.};
    if constexpr(Polynomial::auto_lower_degree)
      EXPECT_EQ(0, p.size());
    else
      EXPECT_EQ(4, p.size());
    p.try_lower_degree();
    EXPECT_EQ(0, p.size());*/
  }

  void test_print() {
    Polynomial p_empty{};
    Polynomial p0{2};
    Polynomial p1{2, 3};
    Polynomial p4{2, 7, 9, 11, 4};

    if constexpr(std::is_same_v<typename Polynomial::coeff_type, double>) {
      EXPECT_PRINT("0", p_empty);
      EXPECT_PRINT("2", p0);
      EXPECT_PRINT("2 + 3*x", p1);
      EXPECT_PRINT("2 + 7*x + 9*x^2 + 11*x^3 + 4*x^4", p4);

    } else {
      EXPECT_PRINT("0", p_empty);
      EXPECT_PRINT("(2,0)", p0);
      EXPECT_PRINT("(2,0) + (3,0)*x", p1);
      EXPECT_PRINT("(2,0) + (7,0)*x + (9,0)*x^2 + (11,0)*x^3 + (4,0)*x^4", p4);
    }
  }

  void test_arithmetics() {
    Polynomial p_empty{};
    Polynomial p1{1, 2, 3, 4};
    Polynomial p2{1, 2, 3, 4, 5, 6};
    Polynomial p3{7, 8, 9, 10};

    // Comparison
    EXPECT_EQ(p1, p1);
    EXPECT_NE(p1, p2);

    // Addition
    EXPECT_EQ(Polynomial({2, 4, 6, 8, 5, 6}), p1 + p2);
    EXPECT_EQ(Polynomial({2, 4, 6, 8, 5, 6}), p2 + p1);
    EXPECT_EQ(Polynomial({8, 10, 12, 14}), p1 + p3);

    // Subtraction
    EXPECT_EQ(Polynomial({0, 0, 0, 0, -5, -6}), p1 - p2);
    EXPECT_EQ(Polynomial({0, 0, 0, 0, 5, 6}), p2 - p1);
    EXPECT_EQ(Polynomial({6, 6, 6, 6}), p3 - p1);

    // Multiplication
    EXPECT_EQ(Polynomial({}), p1 * p_empty);
    EXPECT_EQ(Polynomial({}), p_empty * p2);
    EXPECT_EQ(Polynomial({1, 4, 10, 20, 30, 40, 43, 38, 24}), p1 * p2);
    EXPECT_EQ(Polynomial({}), 2.0 * p_empty);
    EXPECT_EQ(Polynomial({}), p_empty * 2.0);
    EXPECT_EQ(Polynomial({}), 2.0 * p_empty);
    EXPECT_EQ(Polynomial({2, 4, 6, 8}), 2.0 * p1);
    EXPECT_EQ(Polynomial({2, 4, 6, 8}), p1 * 2.0);
  }

  void test_derivative() {
    Polynomial p_empty{};
    Polynomial p0{4};
    Polynomial p5{1, 2, 3, 4, 5, 6};

    EXPECT_EQ(Polynomial({}), derivative(p_empty));
    EXPECT_EQ(Polynomial({}), derivative(p0));
    EXPECT_EQ(Polynomial({2, 6, 12, 20, 30}), derivative(p5));
  }

  void test_antiderivative() {
    Polynomial p_empty{};
    Polynomial p0{4};
    Polynomial p5{1, 2, 3, 4, 5, 6};

    EXPECT_EQ(Polynomial({}), antiderivative(p_empty));
    EXPECT_EQ(Polynomial({0, 4}), antiderivative(p0));
    EXPECT_EQ(Polynomial({0, 1, 1, 1, 1, 1, 1}), antiderivative(p5));
  }
};

using polynomial_types =
    ::testing::Types<polynomial<double, true>, polynomial<double, false>,
                     polynomial<std::complex<double>, true>,
                     polynomial<std::complex<double>, false>>;

#pragma GCC diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
TYPED_TEST_SUITE(polynomials_case, polynomial_types);
TYPED_TEST(polynomials_case, test_empty) { this->test_empty(); }
TYPED_TEST(polynomials_case, test_degree5) { this->test_degree5(); }
TYPED_TEST(polynomials_case, test_try_lower_degree) {
  this->test_try_lower_degree();
}
TYPED_TEST(polynomials_case, test_print) { this->test_print(); }
TYPED_TEST(polynomials_case, test_arithmetics) { this->test_arithmetics(); }
TYPED_TEST(polynomials_case, test_derivative) { this->test_derivative(); }
TYPED_TEST(polynomials_case, test_antiderivative) {
  this->test_antiderivative();
}
