#include <test_tools.hpp>
#include <cmath>

#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#include "spline.hpp"

using std::exp;
using std::sin;
using std::copysign;
using std::abs;
using triqs::arrays::vector;

using namespace som;

// Function to interpolate
double f(double x) { return 2*exp(-0.5*(x-1)*(x-1)) * sin(0.5*x); }

const int N = 20001;
const double xmin = -10, xmax = 10;
const double h = (xmax - xmin)/(N-1);

// Regular mesh
double r_mesh(int i) { return xmin + i*h; }
// Irregular meshes
double ir_mesh1(int i) {
 auto x = r_mesh(i);
 return 0.6*(x + copysign(3 - 3*exp(-2*abs(x)),x));
}
double ir_mesh2(int i) {
 auto x = r_mesh(i);
 return 0.6*(x + copysign(1 - 1*exp(-5*abs(x)),x));
}

TEST(spline,spline) {
 vector<double> x(N), y(N);

 for(int i = 0; i < N; ++i) {
  x(i) = ir_mesh1(i);
  y(i) = f(x(i));
 }

 spline sp(x,y);

 vector<double> ref(N), result(N);
 for(int i = 0; i < N; ++i) {
  double z = r_mesh(i);
  ref(i) = f(z);
  result(i) = sp(z);
 }
 EXPECT_ARRAY_NEAR(result,ref,1e-9);

 for(int i = 0; i < N; ++i) {
  double z = ir_mesh2(i);
  ref(i) = f(z);
  result(i) = sp(z);
 }
 EXPECT_ARRAY_NEAR(result,ref,1e-9);
}

TEST(spline,regular_spline) {
}

MAKE_MAIN;