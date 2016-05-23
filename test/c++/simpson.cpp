#include <triqs/test_tools/arrays.hpp>

#include <cmath>

#include "simpson.hpp"

using namespace som;
using namespace std;

TEST(simpson, exp_sin) {
 auto f = [](double x){ return exp(2*x)*sin(3*x); };

 double ref = exp(6) * (2*sin(9) - 3*cos(9)) / 13;
 ref       -= exp(2) * (2*sin(3) - 3*cos(3)) / 13;
 EXPECT_NEAR(ref, adaptive_simpson(f, 1.0, 3.0, 1e-10), 1e-10);
}

TEST(simpson, sin_inv_x) {
 auto f = [](double x){ return sin(1/x) / (1 + 2*x +x*x); };
 EXPECT_NEAR(0.162985567, adaptive_simpson(f, 0.01, 1.0, 1e-10), 1e-9);
}

MAKE_MAIN;
