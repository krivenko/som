#include <test_tools.hpp>

#include <triqs/arrays.hpp>
#include "triqs/arrays/blas_lapack/gtsv.hpp"

using namespace triqs::arrays;

TEST(gtsv,dgtsv) {

vector<double> DL = {4,3,2,1}; // sub-diagonal elements
vector<double> D = {1,2,3,4,5};// diagonal elements
vector<double> DU = {1,2,3,4}; // super-diagonal elements

vector<double> B1 = {1,2,3,4,5}; // RHS column 1
vector<double> B2 = {6,7,8,9,10}; // RHS column 2
matrix<double> B(5,2);
B(range(),0) = B1;
B(range(),1) = B2;

blas::gtsv(DL,D,DU,B1);
blas::gtsv(DL,D,DU,B2);
blas::gtsv(DL,D,DU,B);

}

TEST(gtsv,cgtsv) {


}

MAKE_MAIN;
