#include <triqs/test_tools/arrays.hpp>

#include "polynomial.hpp"

using som::polynomial;

TEST(polynomial, empty) {
 polynomial p({});

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
